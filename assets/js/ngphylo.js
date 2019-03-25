function init(){
    $('.label-F').addClass('label-success');
    $('.label-P').addClass('label-default');
    $('.label-R').addClass('label-info');
    $('.label-E').addClass('label-danger');
}

$(document).ready(function(){
    init();
});

function toc(){
    var toc="";
    $("h2, h3").each(function(i) {
	var current = $(this);
	current.attr("id", "title" + i);
	toc+="<p class='toc-"+current.prop("tagName")+"'><a id='link" + i + "' href='#title" +
	    i + "' title='" + current.prop("tagName") + "'>" + 
	    current.html() + "</a></p>";
	console.log(current.html());
	current.html('<a class="anchor-title" aria-label="Anchor link for: glyphicons how to use" href="#title'+i+'" class="anchorjs-link ">'+current.html()+'</a>');
    });
    $("#toc").append('<h2>Table of contents</h2>')
    $("#toc").append(toc);
}
