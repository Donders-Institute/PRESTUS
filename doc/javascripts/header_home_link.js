// Make the header title ("PRESTUS - Ultrasound Simulation Toolbox") link to the home page.
// Reuses the logo's href, so the target stays correct under mike version prefixes.
document$.subscribe(function () {
  var logo = document.querySelector(".md-header__button.md-logo");
  var title = document.querySelector(".md-header__title");
  if (!logo || !title || title.dataset.homeLink) return;
  title.dataset.homeLink = "1";
  title.style.cursor = "pointer";
  title.addEventListener("click", function () {
    window.location.href = logo.href;
  });
});
