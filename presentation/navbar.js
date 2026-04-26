// navbar.js — Fortschrittsbalken + Tastaturnavigation
(function() {
  // Fortschrittsbalken
  remark.macros.progress = function() {
    var slideshow = this;
    var el = document.createElement('div');
    el.id = 'progress-bar';
    document.body.appendChild(el);

    slideshow.on('afterShowSlide', function(slide) {
      var total = slideshow.getSlideCount();
      var current = slide.getSlideIndex() + 1;
      el.style.width = (current / total * 100) + '%';
    });
  };

  // Keyboard-Shortcuts anzeigen (? Taste)
  document.addEventListener('keydown', function(e) {
    if (e.key === '?' || e.key === 'h') {
      // Hilfe-Overlay
      var help = document.getElementById('keyboard-help');
      if (help) {
        help.style.display = help.style.display === 'none' ? 'flex' : 'none';
      }
    }
  });
})();