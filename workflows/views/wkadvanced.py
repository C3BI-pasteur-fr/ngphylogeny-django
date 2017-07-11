from django.shortcuts import render, redirect

from galaxy.decorator import connection_galaxy
from tools.models import ToolFlag

WORKFLOW_STATIC_STEPS = [{"step": 0, "category": 'algn', "group": ""},
                         {"step": 1, "category": 'clean', "group": ""},
                         {"step": 2, "category": 'tree', "group": ""},
                         {"step": 3, "category": 'visu', "group": ""},
                        ]

@connection_galaxy
def workflows_advanced_mode_build(request):
    """"""
    for step in WORKFLOW_STATIC_STEPS:
        step['group'] = ToolFlag.objects.get(name=step.get('category'))

    if request.method == 'POST':

        tools = []
        for step in WORKFLOW_STATIC_STEPS:
            tools.append(request.POST.get(step.get('category')))

        request.session['selected_tools'] = tools
        return redirect('workflow_form')

    context = {"workflow": WORKFLOW_STATIC_STEPS}
    return render(request, 'workflows/workflows_advanced.html', context)


