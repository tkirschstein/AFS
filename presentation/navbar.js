/* ============================================================
   navbar.js — Xaringan Presentation Navigation Enhancements
   Fraunhofer IKTS | Agroforestry Biomass Supply Chain Design
   ============================================================

   Features:
     1. Animated top progress bar (IKTS green #179C7D)
     2. Tile-view overlay  (press 'O')
     3. Keyboard shortcut help overlay  (press '?')
     4. Smooth fade-in on slide transition

   Usage in YAML header:
     nature:
       beforeInit: "navbar.js"
   ============================================================ */

'use strict';

(function () {

  /* ── 1. PROGRESS BAR ─────────────────────────────────────── */
  var progressBar = document.createElement('div');
  progressBar.id = 'xar-progress';
  progressBar.style.cssText = [
    'position:fixed',
    'top:0',
    'left:0',
    'height:3px',
    'width:0%',
    'background:#179C7D',
    'z-index:9999',
    'transition:width 0.35s cubic-bezier(0.4,0,0.2,1)',
    'pointer-events:none'
  ].join(';');
  document.body.appendChild(progressBar);

  function updateProgress(slideshow) {
    var current = slideshow.getCurrentSlideIndex() + 1;
    var total   = slideshow.getSlideCount();
    progressBar.style.width = (current / total * 100) + '%';
  }

  /* ── 2. TILE VIEW ────────────────────────────────────────── */
  var tileOverlay = document.createElement('div');
  tileOverlay.id = 'xar-tiles';
  tileOverlay.style.cssText = [
    'display:none',
    'position:fixed',
    'inset:0',
    'background:rgba(0,0,0,0.88)',
    'z-index:9000',
    'overflow-y:auto',
    'padding:2rem'
  ].join(';');

  var tileGrid = document.createElement('div');
  tileGrid.style.cssText = [
    'display:grid',
    'grid-template-columns:repeat(auto-fill,minmax(200px,1fr))',
    'gap:1rem',
    'max-width:1400px',
    'margin:0 auto'
  ].join(';');
  tileOverlay.appendChild(tileGrid);
  document.body.appendChild(tileOverlay);

  var tilesBuilt = false;

  function buildTiles(slideshow) {
    if (tilesBuilt) return;
    tilesBuilt = true;
    var slides = document.querySelectorAll('.remark-slide-container');
    slides.forEach(function (slide, idx) {
      var wrapper = document.createElement('div');
      wrapper.style.cssText = [
        'background:#1e1e1e',
        'border-radius:6px',
        'overflow:hidden',
        'cursor:pointer',
        'border:2px solid transparent',
        'transition:border-color 0.2s'
      ].join(';');
      wrapper.title = 'Folie ' + (idx + 1);

      var num = document.createElement('div');
      num.textContent = idx + 1;
      num.style.cssText = [
        'font-size:0.7rem',
        'color:#aaa',
        'padding:4px 8px',
        'background:#111'
      ].join(';');

      var preview = slide.cloneNode(true);
      preview.style.cssText = [
        'transform:scale(0.25)',
        'transform-origin:top left',
        'width:400%',
        'height:400%',
        'pointer-events:none'
      ].join(';');

      var previewWrap = document.createElement('div');
      previewWrap.style.cssText = 'overflow:hidden;height:120px;position:relative;';
      previewWrap.appendChild(preview);

      wrapper.appendChild(num);
      wrapper.appendChild(previewWrap);

      wrapper.addEventListener('click', function () {
        slideshow.gotoSlide(idx + 1);
        closeTiles();
      });
      wrapper.addEventListener('mouseenter', function () {
        wrapper.style.borderColor = '#179C7D';
      });
      wrapper.addEventListener('mouseleave', function () {
        wrapper.style.borderColor = 'transparent';
      });

      tileGrid.appendChild(wrapper);
    });
  }

  function openTiles(slideshow) {
    buildTiles(slideshow);
    tileOverlay.style.display = 'block';
  }

  function closeTiles() {
    tileOverlay.style.display = 'none';
  }

  tileOverlay.addEventListener('click', function (e) {
    if (e.target === tileOverlay) closeTiles();
  });

  /* ── 3. KEYBOARD HELP OVERLAY ────────────────────────────── */
  var helpOverlay = document.createElement('div');
  helpOverlay.id = 'xar-help';
  helpOverlay.style.cssText = [
    'display:none',
    'position:fixed',
    'inset:0',
    'background:rgba(0,0,0,0.80)',
    'z-index:9100',
    'align-items:center',
    'justify-content:center'
  ].join(';');

  var helpBox = document.createElement('div');
  helpBox.style.cssText = [
    'background:#1e1e1e',
    'color:#e0e0e0',
    'border-radius:10px',
    'padding:2rem 2.5rem',
    'max-width:480px',
    'width:90%',
    'font-family:monospace',
    'font-size:0.9rem',
    'line-height:1.8',
    'border:1px solid #179C7D'
  ].join(';');

  helpBox.innerHTML = [
    '<h3 style="margin:0 0 1rem;color:#179C7D;font-size:1rem;">⌨ Tastenkürzel</h3>',
    '<table style="border-collapse:collapse;width:100%">',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">→ / Space</td><td>Nächste Folie</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">←</td><td>Vorherige Folie</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">O</td><td>Folienübersicht (Tile View)</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">P</td><td>Präsentationsmodus (Notes)</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">F</td><td>Vollbild</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">C</td><td>Klon-Fenster (Presenter)</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">?</td><td>Diese Hilfe</td></tr>',
    '<tr><td style="color:#179C7D;padding-right:1.5rem">Esc</td><td>Schließen / Zurück</td></tr>',
    '</table>',
    '<p style="margin-top:1rem;font-size:0.75rem;color:#888;">Klick außerhalb schließt dieses Fenster.</p>'
  ].join('');

  helpOverlay.appendChild(helpBox);
  document.body.appendChild(helpOverlay);

  function openHelp() {
    helpOverlay.style.display = 'flex';
  }

  function closeHelp() {
    helpOverlay.style.display = 'none';
  }

  helpOverlay.addEventListener('click', function (e) {
    if (e.target === helpOverlay) closeHelp();
  });

  /* ── 4. SLIDESHOW HOOK ───────────────────────────────────── */
  // remark.js fires "ready" after the slideshow is initialised.
  // We hook into it via the global remark macro system.
  var _originalCreate = window.remark && window.remark.create
    ? window.remark.create.bind(window.remark)
    : null;

  if (_originalCreate) {
    window.remark.create = function (options) {
      var slideshow = _originalCreate(options);
      init(slideshow);
      return slideshow;
    };
  }

  // Fallback: wait for the slideshow element to appear in the DOM.
  var _initDone = false;
  function tryInit() {
    if (_initDone) return;
    var ss = window.slideshow || (window.remark && window.remark.slideshow);
    if (ss) { init(ss); _initDone = true; }
  }

  function init(slideshow) {
    // initial progress
    updateProgress(slideshow);

    // update on every slide change
    slideshow.on('showSlide', function () {
      updateProgress(slideshow);
    });

    // keyboard shortcuts
    document.addEventListener('keydown', function (e) {
      var tag = document.activeElement && document.activeElement.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA') return;

      switch (e.key) {
        case 'O': case 'o':
          if (tileOverlay.style.display === 'none') {
            openTiles(slideshow);
          } else {
            closeTiles();
          }
          break;
        case '?':
          if (helpOverlay.style.display === 'none') {
            openHelp();
          } else {
            closeHelp();
          }
          break;
        case 'Escape':
          closeTiles();
          closeHelp();
          break;
      }
    });
  }

  // Try to attach after DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', function () {
      setTimeout(tryInit, 500);
    });
  } else {
    setTimeout(tryInit, 500);
  }

}());
