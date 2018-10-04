import ast
import json
import tempfile

from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy
from django.forms import ValidationError
from django.http import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.views.generic import DetailView, ListView
from bioblend.galaxyclient import ConnectionError
from django.db.utils import OperationalError

from galaxy.decorator import connection_galaxy
from workspace.views import create_history, delete_history
from .forms import ToolForm
from .models import Tool, ToolFieldWhiteList, ToolFlag
from workspace.tasks import monitorworkspace

from Bio import SeqIO, Phylo


class ToolListView(ListView):
    """
        Display the list of executable tool
    """
    CATEGORY = ['algn',
                'clean',
                'tree',
                'visu',
                ]

    def get_queryset(self):
        CATEGORY = ToolFlag.objects.filter(rank=0).values_list(
            'name', flat=True) or self.CATEGORY
        tool_list = Tool.objects.filter(galaxy_server__current=True,
                                        visible=True,
                                        toolflag__name__in=CATEGORY).prefetch_related('toolflag_set')
        return tool_list


class ToolDetailView(DetailView):
    queryset = Tool.objects.filter(galaxy_server__current=True, visible=True)


@connection_galaxy
def tool_exec_view(request, pk, store_output=None):
    """
    :param request:
    :param pk:
    :param store_output:
    :return:
    """

    gi = request.galaxy
    message = ""

    tool_obj = get_object_or_404(Tool, pk=pk)

    toolfield, created = ToolFieldWhiteList.objects.get_or_create(
        tool=tool_obj, context='t')
    tool_visible_field = toolfield.saved_params

    tool_inputs_details = gi.tools.show_tool(
        tool_id=tool_obj.id_galaxy, io_details='true')
    tool_form = ToolForm(tool_params=tool_inputs_details['inputs'],
                         tool_id=pk, whitelist=tool_visible_field,
                         data=request.POST or None,
                         session_files=request.session.get('files')
                         )

    hid = {}
    if request.method == 'POST':
        if tool_form.is_valid():
            try:
                wksph = None 
                tool_inputs = inputs()

                # mapping between form id to galaxy params names
                fields = tool_form.fields_ids_mapping
                exts = tool_form.fields_ext_mapping

                inputs_data = set(tool_form.input_file_ids)

                for key, value in request.POST.items():
                    if not key in inputs_data:
                        tool_inputs.set_param(fields.get(key), value)

                for inputfile in inputs_data:
                    outputs = ""
                    uploaded_file = ""
                    if request.FILES:
                        uploaded_file = request.FILES.get(inputfile, '')
                    if uploaded_file:
                        tmp_file = tempfile.NamedTemporaryFile()
                        for chunk in uploaded_file.chunks():
                            tmp_file.write(chunk)
                        tmp_file.flush()
                        # send file to galaxy after verifying the
                        # allowed extensions
                        type = detect_type(tmp_file.name)
                        if type in ["fasta", "phylip"] and nb_sequences(tmp_file.name, type) <= 3:
                            raise ValueError('Sequence file %s should contain more than 3 sequences for field %s' % (
                                uploaded_file.name, fields.get(inputfile)))
                        if type in exts.get(inputfile, ""):
                            if wksph is None:
                                wksph = create_history(request, name="Analyse with " + tool_obj.name)
                            outputs = gi.tools.upload_file(path=tmp_file.name,
                                                           file_name=uploaded_file.name,
                                                           history_id=wksph.history,
                                                           file_type=type)
                        else:
                            raise ValueError('File format of file %s is not allowed for field %s' % (
                                uploaded_file.name, fields.get(inputfile)))
                    else:
                        # else paste content
                        content = request.POST.get(inputfile)
                        galaxyfile = request.POST.get(
                            "galaxyfile_" + inputfile)
                        if content:
                            tmp_file = tempfile.NamedTemporaryFile()
                            tmp_file.write(content)
                            tmp_file.flush()
                            # send file to galaxy after verifying the
                            # allowed extensions
                            input_fieldname = fields.get(inputfile)
                            type = detect_type(tmp_file.name)
                            if type in exts.get(inputfile, ""):
                                if wksph is None:
                                    wksph = create_history(request, name="Analyse with " + tool_obj.name)
                                outputs = gi.tools.upload_file(path=tmp_file.name,
                                                               file_name=input_fieldname + " pasted_sequence",
                                                               history_id=wksph.history,
                                                               file_type=type)
                            else:
                                raise ValueError(
                                    'This file format is not allowed for field %s' % (input_fieldname))
                        # Else we look at galaxy file ids
                        else:
                            if galaxyfile:
                                file_id = galaxyfile
                                hid[inputfile] = file_id

                    if outputs:
                        file_id = outputs.get('outputs')[0].get('id')
                        hid[inputfile] = file_id

                for h, file_id in hid.items():
                    tool_inputs.set_dataset_param(fields.get(h), file_id)

                if wksph is None:
                    wksph = create_history(request, name="Analyse with " + tool_obj.name)
                    
                tool_outputs = gi.tools.run_tool(history_id=wksph.history,
                                                 tool_id=tool_obj.id_galaxy,
                                                 tool_inputs=tool_inputs)
                # Start monitoring (for sending emails)
                monitorworkspace.delay(wksph.history)
                wksph.monitored = True
                wksph.save()

                if store_output:
                    request.session['output'] = tool_outputs

                return redirect(reverse_lazy("history_detail", kwargs={'history_id': wksph.history}))

            except ValueError as ve:
                message = str(ve)
            except ConnectionError as ce:
                message = str(ce)
            except NameError as ne:
                message = str(ne)
            except OperationalError as oe:
                message = str(oe)
            except AttributeError as ae:
                message = str(ae)
            except TypeError as te:
                message = str(te)
            except Exception as e:
                raw_message = ast.literal_eval(e.body)
                reverse_dict_field = {
                    v: k for k, v in tool_form.fields_ids_mapping.items()}
                err_data = raw_message.get('err_data')

                if err_data:
                    for field, err_msg in err_data.items():
                        tool_form.add_error(reverse_dict_field.get(field),
                                            ValidationError(err_msg, code='invalid'))
                else:
                    message = raw_message

                delete_history(request, wksph.history)

    context = {"toolform": tool_form,
               "tool": tool_obj,
               "message": message}

    return render(request, 'tools/tool_form.html', context)


@connection_galaxy
def get_tool_name(request):
    """
    Find the name of tool from an id_tool (used by AJAX)
    """
    context = dict()

    if request.POST:
        gi = request.galaxy
        toolid = request.POST.get('tool_id')

        if toolid:
            tool = gi.tools.get_tools(tool_id=toolid)[0]
            context.update({'tool_id': toolid, 'name': tool.get('name')})

    return HttpResponse(json.dumps(context), content_type='application/json')


def detect_type(filename):
    """
    :param filename: File to read and detect the format
    :return: detected type, in [fasta, phylip, phylip-relaxed, newick, N/A]

    Tests formats using biopython SeqIO or Phylo
    """
    # Check Fasta Format
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "fasta"):
            nbseq += 1
        if nbseq > 0:
            return "fasta"
    except Exception:
        pass

    # Check phylip strict
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "phylip"):
            nbseq += 1
        if nbseq > 0:
            return "phylip"
    except Exception:
        pass

    # Check phylip relaxed
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "phylip-relaxed"):
            nbseq += 1
        if nbseq > 0:
            return "phylip"
    except Exception:
        pass

    # Check Newick
    try:
        nbtrees = 0
        trees = Phylo.parse(filename, 'newick')
        for t in trees:
            nbtrees += 1
        if nbtrees > 0:
            return "nhx"
    except Exception as e:
        pass

    return "txt"


def nb_sequences(filename, format):
    nbseq = 0
    if format == 'fasta':
        try:
            for r in SeqIO.parse(filename, "fasta"):
                nbseq += 1
        except Exception:
            pass
    elif format == 'phylip':
        try:
            for r in SeqIO.parse(filename, "phylip"):
                nbseq += 1
        except Exception:
            pass

    return nbseq
