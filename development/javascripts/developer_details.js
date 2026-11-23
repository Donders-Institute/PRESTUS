// Header switch "Developer details": shows or hides `developer` boxes
// (??? developer / !!! developer) and the Developers nav section on every page.
// Off by default; the choice is kept in localStorage. The hiding itself is in
// stylesheets/extra.css (html[data-dev="on"]), so pages do not flash on load.
(function () {
  var KEY = "prestus-developer-details";

  function load() {
    try { return localStorage.getItem(KEY) === "on"; } catch (e) { return false; }
  }
  function save(on) {
    try { localStorage.setItem(KEY, on ? "on" : "off"); } catch (e) { /* private mode */ }
  }
  function apply(on) {
    document.documentElement.dataset.dev = on ? "on" : "off";
    var btn = document.querySelector(".prestus-dev-toggle");
    if (btn) {
      btn.setAttribute("aria-pressed", on ? "true" : "false");
      btn.title = on ? "Hide developer details" : "Show developer details";
    }
    var note = document.querySelector(".prestus-dev-notice");
    if (note && on) note.remove();
  }

  apply(load());

  document$.subscribe(function () {
    // Mark the Developers nav section so the stylesheet can hide it.
    document.querySelectorAll(".md-nav--primary > .md-nav__list > .md-nav__item").forEach(function (item) {
      var label = item.querySelector(":scope > .md-nav__link, :scope > label");
      if (label && /Developers/.test(label.textContent)) item.classList.add("prestus-dev-nav");
    });

    // Headings marked { .developer }: hide their table-of-contents entries too.
    document.querySelectorAll(".md-typeset h2.developer[id], .md-typeset h3.developer[id]").forEach(function (h) {
      document.querySelectorAll('.md-nav--secondary a[href="#' + h.id + '"]').forEach(function (a) {
        a.parentElement.classList.add("prestus-dev-nav");
      });
    });

    var header = document.querySelector(".md-header__inner");
    if (header && !header.querySelector(".prestus-dev-toggle")) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className = "md-header__button md-icon prestus-dev-toggle";
      btn.setAttribute("aria-label", "Developer details");
      btn.innerHTML =
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
        '<path d="M14.6 16.6 19.2 12l-4.6-4.6L16 6l6 6-6 6zm-5.2 0L4.8 12l4.6-4.6L8 6l-6 6 6 6z"/></svg>' +
        '<span class="prestus-dev-toggle__label">Developer</span>';
      btn.addEventListener("click", function () {
        var on = document.documentElement.dataset.dev !== "on";
        save(on);
        apply(on);
      });
      var anchor = header.querySelector(".md-search") || header.lastElementChild;
      header.insertBefore(btn, anchor);
    }
    apply(load());

    // A link or search result pointing into a hidden developer box: offer to show it.
    var id = decodeURIComponent(location.hash.slice(1));
    var target = id && document.getElementById(id);
    if (target && target.closest(".developer") && document.documentElement.dataset.dev !== "on"
        && !document.querySelector(".prestus-dev-notice")) {
      var notice = document.createElement("div");
      notice.className = "admonition prestus-dev-notice";   // not .developer: that is hidden
      notice.innerHTML = '<p class="admonition-title">This section is developer detail</p>' +
        '<p><button type="button" class="md-button">Show developer details</button></p>';
      notice.querySelector("button").addEventListener("click", function () {
        save(true);
        apply(true);
        target.scrollIntoView();
      });
      var content = document.querySelector(".md-content__inner");
      if (content) content.insertBefore(notice, content.firstChild);
    }
  });
})();
