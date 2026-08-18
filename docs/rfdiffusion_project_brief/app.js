(function () {
  "use strict";

  const params = new URLSearchParams(window.location.search);
  const requestedView = params.get("view") === "scroll" ? "scroll" : "slide";

  document.querySelectorAll("[data-view]").forEach((link) => {
    const isActive = link.dataset.view === (requestedView === "scroll" ? "scroll" : "slides");
    link.classList.toggle("active", isActive);
    if (isActive) link.setAttribute("aria-current", "page");
  });

  Reveal.initialize({
    hash: true,
    view: requestedView,
    scrollProgress: true,
    scrollSnap: "proximity",
    scrollLayout: "full",
    controls: true,
    controlsTutorial: true,
    progress: true,
    slideNumber: "c/t",
    showSlideNumber: "all",
    center: false,
    width: 1440,
    height: 900,
    margin: 0.035,
    minScale: 0.2,
    maxScale: 1.5,
    transition: "fade",
    backgroundTransition: "fade",
    pdfSeparateFragments: false,
    plugins: [RevealNotes]
  });
})();
