var branch_length_accept = true;
var zoom_value = 1;
var presto = {};

presto.zoom_value = 1;
presto.branch_length_accept = true;

load = function(newick_content) {

    function  mySVGCheckFn()
    {
        var  bbox = $('svg')[0].getBBox();

        $('svg').attr("width", bbox.x + bbox.width + 70);
        $('svg').attr("height",bbox.y + bbox.height + 100);
    }

    $("#newick_file").on("change", function (e) {
        var files = e.target.files; // FileList object

        if (files.length == 1) {
            var f = files[0];
            var reader = new FileReader();

            reader.onload = (function (theFile) {
                return function (e) {
                    var res = d3_phylotree_newick_parser(e.target.result);

                    if (res["json"]) {
                        if (!("children" in res["json"])) {
                            res["error"] = "Empty tree";
                        }
                    }

                    var warning_div = d3.select("#main_display").insert("div", ":first-child");
                    if (res["error"]) {
                        warning_div.attr("class", "alert alert-danger alert-dismissable")
                            .html("<strong>Newick parser error for file " + f.name + ": </strong> In file " + res["error"]);

                    } else {
                        default_tree_settings();
                        tree(res).svg(treeSvg).layout();
                        warning_div.attr("class", "alert alert-success alert-dismissable")
                            .html("Loaded a tree from  file <strong>" + f.name + ": </strong>");
                    }
                    warning_div.append("button")
                        .attr("type", "button")
                        .attr("class", "close")
                        .attr("data-dismiss", "alert")
                        .attr("aria-hidden", "true")
                        .html("&times;");
                };
            })(f);

            reader.readAsText(f);
        }
    });

    $("#mp_label").on("click", function (e) {
        tree.max_parsimony(true);
    });

    $("[data-direction]").on("click", function (e) {
        var which_function = $(this).data("direction") == 'vertical' ? tree.spacing_x : tree.spacing_y;
        which_function(which_function() + (+$(this).data("amount"))).update();
        mySVGCheckFn();
    });


    $(".phylotree-layout-mode").on("change", function (e) {
        if ($(this).is(':checked')) {
            if (tree.radial() != ($(this).data("mode") == "radial")) {
                tree.radial(!tree.radial()).placenodes().update();
            }
            if ($(this).data("mode") == "straight") {
                tree.options({'branches': 'straight'}, true);
            }
            else if ($(this).data("mode") == "step") {
                tree.options({'branches': 'step'}, true);
            }
        }
    });

    $(".phylotree_branch_length").on("change", function (e) { // HERE
        if ($(this).is(':checked')) {
            if (branch_length_accept == true) {
                branch_length_accept = false;
                if (typeof new_json != 'undefined') {
                }
                tree.placenodes().update();
            }
            else {
                branch_length_accept = true;
                if (typeof new_json != 'undefined') {
                }
                tree.placenodes().update();
            }
        }
    });

    $("#align-toggler").on("click", function (e) { // HERE : a changer ?
        if ($('#align-toggler').is(":checked")) {
            tree.options({'align-tips': 'right'}, true);
        }
        else {
            tree.options({'align-tips': false}, true);
        }
    });

    function sort_nodes(asc) {
        tree.traverse_and_compute(function (n) {
            var d = 1;
            if (n.children && n.children.length) {
                d += d3.max(n.children, function (d) {
                    return d["count_depth"];
                });
            }
            n["count_depth"] = d;
        });
        tree.resort_children(function (a, b) {
            return (a["count_depth"] - b["count_depth"]) * (asc ? 1 : -1);
        });
    }

    $("#sort_original").on("click", function (e) {
        tree.resort_children(function (a, b) {
            return a["original_child_order"] - b["original_child_order"];
        });
    });

    $("#sort_ascending").on("click", function (e) {
        sort_nodes(true);
    });

    $("#sort_descending").on("click", function (e) {
        sort_nodes(false);
    });

    $("#branch_filter").on("input propertychange", function (e) {
        var filter_value = $(this).val();

        var rx = new RegExp(filter_value, "i");

        tree.modify_selection(function (n) {
            return filter_value.length && (tree.branch_name()(n.target).search(rx)) != -1;
        }, "tag");

    });

    $("#validate_newick").on("click", function (e) {
        var res = d3_phylotree_newick_parser($('textarea[id$="nwk_spec"]').val(), true);
        if (res["error"] || !res["json"]) {
            var warning_div = d3.select("#newick_body").selectAll("div  .alert-danger").data([res["error"]])
            warning_div.enter().append("div");
            warning_div.html(function (d) {
                return d;
            }).attr("class", "alert-danger");

        } else {
            default_tree_settings();
            tree(res).svg(treeSvg).layout();
            $('#newick_modal').modal('hide');
        }
    });

    function default_tree_settings() {
        tree.branch_length(null);
        tree.branch_name(null);
        tree.node_span('equal');
        tree.options({'draw-size-bubbles': false}, false);
        tree.style_nodes();
        tree.style_edges();
    }

    function update_controls() {
        $("[data-mode='" + (tree.radial() ? 'radial' : 'linear') + "']").click();
        $("[data-align='" + (tree.align_tips() ? 'right' : 'left') + "']").click();
    }

    var display_bootstrap = false;
    var display_LB = false;

    $("#display_bootstrap").on("click", function (e) { // HERE : fonctionnalité d'affichage bootstrap / LB CHANGER DE PLACE ?
        if ($('#display_bootstrap').is(":checked")) {
            display_bootstrap = true;
        }
        else {
            display_bootstrap = false;
        }
        tree.style_nodes(node_annotater);
        tree.update();
    });

    $("#display_LB").on("click", function (e) {
        if ($('#display_LB').is(":checked")) {
            display_LB = true;
        }
        else {
            display_LB = false;
        }
        tree.style_nodes(node_annotater);
        tree.update();
    });

    function node_annotater(element, data) { // HERE EN TRAVAUX
        if ("bootstrap_values" in data && data.bootstrap_values) {
            var label = element.selectAll(".bootstrap_values");
            if (label.empty() && display_bootstrap == true) {
                if (data.bootstrap_values.length >= 4) {
                    element.append("text").classed("bootstrap_values", true).text(data.bootstrap_values.substring(0, 4)).attr("dx", ".3em").attr("text-anchor", "start").attr("alignment-baseline", "middle");

                }
                else {
                    element.append("text").classed("bootstrap_values", true).text(data.bootstrap_values).attr("dx", ".3em").attr("text-anchor", "start").attr("alignment-baseline", "middle");
                }
            }
            else if (!label.empty() && display_bootstrap == false) {
                label.remove();
            }
        }
        if ("attribute" in data && data.attribute) {
            var branch_length = element.selectAll(".attribute");
            if (branch_length.empty() && display_LB == true) {
                if (data.attribute.length >= 4) {
                    element.append("text").classed("attribute", true).text(data.attribute.substring(0, 4)).attr("dx", "-.6em").attr("text-anchor", "end").attr("alignment-baseline", "middle");
                }
                else {
                    element.append("text").classed("attribute", true).text(data.attribute).attr("dx", "-.6em").attr("text-anchor", "end").attr("alignment-baseline", "middle");

                }
            }
            else if (!branch_length.empty() && display_LB == false) {
                branch_length.remove();
            }
        }
    }
    var width = 790,
        height = 850, //$(container_id).height()
        selection_set = ['Foreground'],
        current_selection_name = $("#selection_name_box").val(),
        current_selection_id = 0,
        max_selections = 10;
    color_scheme = d3.scale.category10(),
        selection_menu_element_action = "phylotree_menu_element_action";

    var tree = d3.layout.phylotree("body")
        .size([height, width])
        .separation(function (a, b) {
            return 0;
        })
        .count_handler(function (count) {
                $("#selected_branch_counter").text(function (d) {
                    return count[current_selection_name];
                });
                $("#selected_filtered_counter").text(count.tag);
            }
        );

    var text_nwk = d3_phylotree_newick_parser(newick_content, true);

    var container_id = '#tree_container';

    var zoom = d3.behavior.zoom() // fonction de zoom
        .size([100])
        .scaleExtent([0.5, 3])
        .scale(1.2)
        .on("zoom", zoomed);

    var enable_zoom = false;

    $("#enable_zoom").on("click", function (e) {
        if (enable_zoom == false) {
            enable_zoom = true;
        }
        else {
            enable_zoom = false;
        }
    });

    $("#zoom_plus").on("click", function (e) {
        if ((zoom_value * 1.2) <= 3) {
            if ((zoom_value * 1.2) <= 1) {
                zoom_value *= 1.2;
                $('.phylotree-container').attr("transform", "translate(" + presto.translate_value_tree + ")scale(" + zoom_value + ")");
                $('.tree-scale-bar').attr("transform", "translate(" + presto.translate_value_scale + ")scale(" + zoom_value + ")");
                mySVGCheckFn();
            }
            else {
                zoom_value *= 1.2;
                $('.phylotree-container').attr("transform", "translate(" + (presto.translate_value_tree[0] + (10*zoom_value)) + "," + (presto.translate_value_tree[1] + (10*zoom_value)) + ")scale(" + zoom_value + ")");
                $('.tree-scale-bar').attr("transform", "translate(" + (presto.translate_value_scale[0] + (10*zoom_value)) + "," + (presto.translate_value_scale[1] + (10*zoom_value)) +  ")scale(" + zoom_value + ")");
                mySVGCheckFn();
            }
        }
        else {
            zoom_value = 3;
            $('.phylotree-container').attr("transform", "translate(" + (presto.translate_value_tree[0] + (10*zoom_value)) + "," + (presto.translate_value_tree[1] + (10*zoom_value)) + ")scale(" + zoom_value + ")");
            $('.tree-scale-bar').attr("transform", "translate(" + (presto.translate_value_scale[0] + (10*zoom_value)) + "," + (presto.translate_value_scale[1] + (10*zoom_value)) +  ")scale(" + zoom_value + ")");
            mySVGCheckFn();
        }

    });

    $("#zoom_minus").on("click", function (e) {
        if ((zoom_value /1.2) >= 0.5){
            if (((zoom_value / 1.2) >= 1)) {
                zoom_value /= 1.2;
                $('.phylotree-container').attr("transform", "translate(" + (presto.translate_value_tree[0] + (10*zoom_value)) + "," + (presto.translate_value_tree[1] + (10*zoom_value)) + ")scale(" + zoom_value + ")");
                $('.tree-scale-bar').attr("transform", "translate(" + (presto.translate_value_scale[0] + (10*zoom_value)) + "," + (presto.translate_value_scale[1] + (10*zoom_value)) +  ")scale(" + zoom_value + ")");
                mySVGCheckFn();
            }
            else {
                zoom_value /= 1.2;
                $('.phylotree-container').attr("transform", "translate(" + presto.translate_value_tree + ")scale(" + zoom_value + ")");
                $('.tree-scale-bar').attr("transform", "translate(" + presto.translate_value_scale + ")scale(" + zoom_value + ")");
                mySVGCheckFn();
            }
        }
        else {
            zoom_value = 0.5;
            $('.phylotree-container').attr("transform", "translate(" + presto.translate_value_tree + ")scale(" + zoom_value + ")");
            $('.tree-scale-bar').attr("transform", "translate(" + presto.translate_value_scale + ")scale(" + zoom_value + ")");
            mySVGCheckFn();
        }
    });



    function zoomed() {
        if (enable_zoom == true) {
            zoom_value = d3.event.scale;
            if (zoom_value <= 1) {

                $('.phylotree-container').attr("transform", "translate(" + presto.translate_value_tree + ")scale(" + zoom_value + ")");
                $('.tree-scale-bar').attr("transform", "translate(" + presto.translate_value_scale + ")scale(" + zoom_value + ")");
                mySVGCheckFn();
            }
            else {
                $('.phylotree-container').attr("transform", "translate(" + (presto.translate_value_tree[0] + (10*zoom_value)) + "," + (presto.translate_value_tree[1] + (10*zoom_value)) + ")scale(" + zoom_value + ")");
                $('.tree-scale-bar').attr("transform", "translate(" + (presto.translate_value_scale[0] + (10*zoom_value)) + "," + (presto.translate_value_scale[1] + (10*zoom_value)) +  ")scale(" + zoom_value + ")");
                mySVGCheckFn();
            }
        }

    }

    var treeSvg = d3.select(container_id).append("svg")
        .call(zoom);

    $(document).ready(function () {
        default_tree_settings();
        tree(text_nwk).svg(treeSvg).layout();
    });

    $("#download-SVG").on("click", function (e) { // HERE : choix du format
        var format = "svg";
        var content = $("#tree_container").html();
        var formulaire = $("<form id='formulaire_download' action='exportSVG.php' method='post'> <input type='hidden' name= 'content' value = '" + content + "'> <input type='hidden' name= 'format' value = '" + format + "' ><input type='hidden' name = 'file_name' value = 'file' ></form>");
        $("body").append(formulaire);
        var form_get = document.getElementById('formulaire_download');
        form_get.submit();
        form_get.parentNode.removeChild(form_get);
    });

    $("#download-PDF").on("click", function (e) {
        var format = "pdf";
        var content = $("#tree_container").html();
        var formulaire = $("<form id='formulaire_download' action='exportSVG.php' method='post'> <input type='hidden' name= 'content' value = '" + content + "'> <input type='hidden' name= 'format' value = '" + format + "' ><input type='hidden' name = 'file_name' value = 'file' ></form>");
        $("body").append(formulaire);
        var form_get = document.getElementById('formulaire_download');
        form_get.submit();
        form_get.parentNode.removeChild(form_get);
    });

    $("#download-PNG").on("click", function (e) {
        var format = "png";
        var content = $("#tree_container").html();
        var formulaire = $("<form id='formulaire_download' action='exportSVG.php' method='post'> <input type='hidden' name= 'content' value = '" + content + "'> <input type='hidden' name= 'format' value = '" + format + "' ><input type='hidden' name = 'file_name' value = 'file' ></form>");
        $("body").append(formulaire);
        var form_get = document.getElementById('formulaire_download');
        form_get.submit();
        form_get.parentNode.removeChild(form_get);
    });

    $("#download-JPEG").on("click", function (e) {
        var format = "jpeg";
        var content = $("#tree_container").html();
        var formulaire = $("<form id='formulaire_download' action='exportSVG.php' method='post'> <input type='hidden' name= 'content' value = '" + content + "'> <input type='hidden' name= 'format' value = '" + format + "' ><input type='hidden' name = 'file_name' value = 'file' ></form>");
        $("body").append(formulaire);
        var form_get = document.getElementById('formulaire_download');
        form_get.submit();
        form_get.parentNode.removeChild(form_get);
    });

    $("#download-NWK").on("click", function (e) {
        var format = "nwk";
        var content = newick_content;
        var formulaire = $("<form id='formulaire_download' action='exportSVG.php' method='post'> <input type='hidden' name= 'content' value = '" + content + "'> <input type='hidden' name= 'format' value = '" + format + "' ><input type='hidden' name = 'file_name' value = 'file' ></form>");
        $("body").append(formulaire);
        var form_get = document.getElementById('formulaire_download');
        form_get.submit();
        form_get.parentNode.removeChild(form_get);
    });
}
