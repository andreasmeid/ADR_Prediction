(function () {
var script = document.createElement("script");
script.type = "text/javascript";
script.src  = "libs/mathjax-local/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
document.getElementsByTagName("head")[0].appendChild(script);
})();

remark.macros.scale = function (percentage) {
  var url = this;
  return '<img src="' + url + '" style="width: ' + percentage + '" />';
};
