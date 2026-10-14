/* SOLARIS · Privacidad: inventario de lo que se guarda en este navegador, preferencia de contenido externo (fotografías de
   Wikimedia Commons, desactivada por defecto) y panel «Gestionar privacidad y cookies». Lo comparten el simulador y las
   páginas legales. Nada de esto se envía a ningún servidor: SOLARIS es un sitio estático sin servidor propio. */
(function () {
'use strict';
var KEY = 'solaris.privacidad', VER = 1, MAX_HIST = 12;
// versiones vigentes de los documentos legales (cámbialas aquí y en cada página cuando se actualice su texto)
var DOCS = {
  terminos: { v: '1.1', titulo: 'Términos y condiciones', ruta: 'terminos-y-condiciones/' },
  privacidad: { v: '1.1', titulo: 'Política de privacidad', ruta: 'politica-de-privacidad/' },
  cookies: { v: '1.1', titulo: 'Política de cookies', ruta: 'politica-de-cookies/' },
};
// raíz del sitio a partir de la ubicación de este script (sirve igual desde el simulador y desde las subpáginas)
var ROOT = (function () { try { var s = document.currentScript && document.currentScript.src; return s ? new URL('..', s).href : './'; } catch (e) { return './'; } })();
// almacenamiento local que SOLARIS puede crear; ninguna clave contiene datos que identifiquen a una persona
var STORE = [
  { key: 'solaris.settings', nombre: 'Ajustes de SOLARIS', fin: 'Tus ajustes de gráficos, interfaz, controles, audio, simulación y accesibilidad, la última fecha simulada y si ya viste la bienvenida, para que se conserven entre visitas.' },
  { key: 'solaris.vuelo', nombre: 'Progreso del modo de vuelo', fin: 'Bitácora de objetos escaneados, logros, mejores tiempos de los retos, estadísticas de vuelo, si completaste la escuela de vuelo y si la voz de a bordo está activada.' },
  { key: KEY, nombre: 'Preferencias de privacidad', fin: 'Si permitiste cargar fotografías de Wikimedia Commons y el historial local de esa decisión.' },
  { key: 'solaris.music', nombre: 'Música (formato antiguo)', fin: 'Preferencia de música de versiones anteriores. Ya no se escribe: solo se lee una vez para trasladarla a los ajustes.', antigua: true },
  { key: 'solaris.audio', nombre: 'Audio (formato antiguo)', fin: 'Volúmenes de versiones anteriores. Ya no se escribe: solo se lee una vez para trasladarlos a los ajustes.', antigua: true },
];
var PURPOSES = {
  wikimedia: { nombre: 'Fotografías reales de Wikimedia Commons', tipo: 'Contenido externo opcional' },
};
var METODOS = { panel: 'Panel de privacidad', ficha: 'Botón en la ficha de un objeto', comparador: 'Botón del comparador', visor: 'Visor de fotografías', borrado: 'Borrado de datos' };

function ls() { try { return window.localStorage; } catch (e) { return null; } }
function esc(s) { return String(s).replace(/[&<>"']/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]; }); }
function blank() { return { v: VER, ext: { wikimedia: false }, hist: [] }; }
function read() {
  var s = ls(), d = blank();
  if (!s) return d;
  try {
    var o = JSON.parse(s.getItem(KEY) || 'null');
    if (o && o.v === VER) {
      if (o.ext && typeof o.ext === 'object') for (var k in d.ext) d.ext[k] = o.ext[k] === true;
      if (Array.isArray(o.hist)) d.hist = o.hist.filter(function (h) { return h && typeof h === 'object' && PURPOSES[h.fin]; }).slice(-MAX_HIST);
    }
  } catch (e) { }
  return d;
}
function write(d) { var s = ls(); if (!s) return false; try { s.setItem(KEY, JSON.stringify(d)); return true; } catch (e) { return false; } }
function emit() { var d = read(); try { window.dispatchEvent(new CustomEvent('solaris:privacy', { detail: { wikimedia: d.ext.wikimedia } })); } catch (e) { } }
function now() { return new Date().toISOString().replace(/\.\d{3}Z$/, 'Z'); }
function when(iso) { try { return new Date(iso).toLocaleString('es-MX', { dateStyle: 'medium', timeStyle: 'short' }); } catch (e) { return iso; } }
function bytes(k) { var s = ls(); if (!s) return 0; try { var v = s.getItem(k); return v == null ? -1 : (k.length + v.length) * 2; } catch (e) { return 0; } }
function size(n) { return n < 1024 ? n + ' B' : (n / 1024).toFixed(n < 10240 ? 1 : 0).replace('.', ',') + ' KB'; }

var P = {
  KEY: KEY, DOCS: DOCS, ROOT: ROOT, STORE: STORE,
  gpc: function () { return navigator.globalPrivacyControl === true; },
  allows: function (fin) { return read().ext[fin] === true; },
  // decisión del usuario: activa o retira un permiso y deja constancia solo en este navegador
  set: function (fin, ok, metodo) {
    if (!PURPOSES[fin]) return false;
    var d = read(), t = now();
    if (d.ext[fin] === !!ok) return true;
    d.ext[fin] = !!ok;
    if (!ok) d.hist.forEach(function (h) { if (h.fin === fin && h.ok && !h.retirado) h.retirado = t; });
    d.hist.push({ id: 'd' + Date.now().toString(36), fin: fin, ok: !!ok, doc: 'Política de privacidad ' + DOCS.privacidad.v, t: t, metodo: METODOS[metodo] ? metodo : 'panel', retirado: null });
    d.hist = d.hist.slice(-MAX_HIST);
    var saved = write(d);
    emit(); P.refresh();
    return saved;
  },
  // borra claves de SOLARIS de este navegador; en el simulador avisa antes para que no se vuelvan a escribir y recarga
  wipe: function (keys) {
    var s = ls(); if (!s) return;
    var app = !!document.getElementById('app');
    try { window.dispatchEvent(new CustomEvent('solaris:privacy-wipe', { detail: { keys: keys.slice() } })); } catch (e) { }
    keys.forEach(function (k) { try { s.removeItem(k); } catch (e) { } });
    if (keys.indexOf(KEY) >= 0) emit();
    if (app && keys.some(function (k) { return k !== KEY; })) { P.say('Datos borrados. Recargando SOLARIS con los valores predeterminados…'); setTimeout(function () { location.reload(); }, 700); return; }
    P.say(keys.length > 1 ? 'Se borraron los datos de SOLARIS de este navegador.' : 'Dato borrado de este navegador.');
    P.refresh();
  },
  url: function (doc, hash) { return ROOT + DOCS[doc].ruta + (hash ? '#' + hash : ''); },
  linksHTML: function () {
    return '<a href="' + P.url('terminos') + '" target="_blank" rel="noopener">Términos y condiciones</a>'
      + '<a href="' + P.url('privacidad') + '" target="_blank" rel="noopener">Política de privacidad</a>'
      + '<a href="' + P.url('cookies') + '" target="_blank" rel="noopener">Política de cookies</a>';
  },
  // ---------- panel ----------
  el: null, from: null,
  build: function () {
    if (P.el) return P.el;
    var w = document.createElement('div');
    w.id = 'pv-panel'; w.className = 'pv-wrap'; w.hidden = true;
    w.setAttribute('role', 'dialog'); w.setAttribute('aria-modal', 'true'); w.setAttribute('aria-labelledby', 'pv-title');
    w.innerHTML = '<div class="pv-card">'
      + '<header class="pv-head"><div><span class="pv-kicker">SOLARIS · Privacidad</span><h2 id="pv-title">Gestionar privacidad y cookies</h2></div>'
      + '<button type="button" class="pv-x" data-pv-close aria-label="Cerrar">&times;</button></header>'
      + '<div class="pv-body" tabindex="-1">'
      + '<p class="pv-lead">SOLARIS <b>no usa cookies</b>, no tiene cuentas ni formularios y no incluye analítica, publicidad ni herramientas de rastreo. Tampoco tiene un servidor propio que reciba tus datos.</p>'
      + '<section aria-labelledby="pv-h1"><h3 id="pv-h1">Datos guardados en este navegador</h3>'
      + '<p class="pv-note">Solo se crean cuando usas la función correspondiente. Se quedan en tu dispositivo, no caducan solos, no se envían a nadie y no sirven para identificarte. Son necesarios para darte el servicio que pides, por eso no se solicita consentimiento; aun así puedes borrarlos cuando quieras.</p>'
      + '<ul class="pv-store" id="pv-store"></ul>'
      + '<button type="button" class="pv-btn pv-danger" data-pv-wipe="*">Borrar todos los datos de SOLARIS en este navegador</button></section>'
      + '<section aria-labelledby="pv-h2"><h3 id="pv-h2">Contenido externo opcional</h3>'
      + '<div class="pv-opt"><label class="pv-sw"><input type="checkbox" role="switch" id="pv-wiki" aria-describedby="pv-wiki-d"><span class="pv-track" aria-hidden="true"></span><b>Fotografías reales de Wikimedia Commons</b></label>'
      + '<p id="pv-wiki-d" class="pv-note">Muestra fotografías reales junto a los modelos en las fichas, el comparador y el modo de vuelo. Al activarlo, tu navegador descarga las imágenes directamente de servidores de la Fundación Wikimedia (Estados Unidos), que recibe tu dirección IP y los datos técnicos habituales de una conexión web conforme a su <a href="https://foundation.wikimedia.org/wiki/Policy:Privacy_policy" target="_blank" rel="noopener">política de privacidad</a>. SOLARIS no recibe esos datos. Está desactivado por defecto y puedes retirarlo en cualquier momento: las fotos dejan de cargarse al instante.</p>'
      + '<p class="pv-state" id="pv-wiki-s" aria-live="polite"></p></div>'
      + '<p class="pv-gpc" id="pv-gpc" hidden>Tu navegador envía la señal <b>Global Privacy Control</b>. SOLARIS no vende ni comparte datos personales, así que no hay nada que excluir; el contenido externo sigue desactivado salvo que tú lo actives.</p></section>'
      + '<section aria-labelledby="pv-h3"><h3 id="pv-h3">Lo que no depende de este panel</h3><ul class="pv-list">'
      + '<li><b>Alojamiento.</b> SOLARIS se publica en GitHub Pages. Para entregarte la página, GitHub recibe tu dirección IP y la registra por motivos de seguridad. Es necesario para que el sitio funcione y lo gestiona GitHub conforme a su <a href="https://docs.github.com/es/site-policy/privacy-policies/github-general-privacy-statement" target="_blank" rel="noopener">declaración de privacidad</a>.</li>'
      + '<li><b>Tipografías, texturas, modelos y música</b> se sirven desde el propio sitio: no se conectan a servicios de terceros.</li>'
      + '<li><b>Voz de a bordo.</b> Usa solo voces instaladas en tu dispositivo; las voces en línea del navegador no se utilizan.</li>'
      + '<li><b>Enlaces externos</b> (créditos, licencias, fuentes científicas): solo te conectan con otro sitio si los abres, y ese sitio aplica sus propias políticas.</li>'
      + '</ul></section>'
      + '<section aria-labelledby="pv-h4"><h3 id="pv-h4">Historial de tus decisiones</h3><ol class="pv-hist" id="pv-hist"></ol>'
      + '<p class="pv-note">Se guarda solo en este navegador para que puedas revisar lo que decidiste. No se envía a nadie ni es un registro verificable: desaparece si borras los datos del navegador.</p></section>'
      + '<p class="pv-msg" id="pv-msg" role="status" aria-live="polite"></p>'
      + '</div>'
      + '<footer class="pv-foot"><nav class="pv-docs" aria-label="Documentos legales">' + P.linksHTML() + '</nav><button type="button" class="pv-btn pv-ok" data-pv-close>Listo</button></footer>'
      + '</div>';
    document.body.appendChild(w);
    P.el = w;
    w.addEventListener('click', function (e) {
      if (e.target === w || e.target.closest('[data-pv-close]')) { P.close(); return; }
      var b = e.target.closest('[data-pv-wipe]'); if (!b) return;
      if (b.dataset.armed !== '1') {           // segunda pulsación para confirmar: evita borrados accidentales
        b.dataset.armed = '1'; b.dataset.txt = b.textContent;
        b.textContent = document.getElementById('app') && b.dataset.pvWipe !== KEY ? 'Confirmar: borrar y recargar' : 'Confirmar borrado';
        clearTimeout(b._t); b._t = setTimeout(function () { b.dataset.armed = ''; b.textContent = b.dataset.txt; }, 5000);
        return;
      }
      clearTimeout(b._t);
      P.wipe(b.dataset.pvWipe === '*' ? STORE.map(function (s) { return s.key; }) : [b.dataset.pvWipe]);
    });
    w.addEventListener('change', function (e) { if (e.target.id === 'pv-wiki') P.set('wikimedia', e.target.checked, 'panel'); });
    return w;
  },
  refresh: function () {
    if (!P.el) return;
    var d = read(), st = document.getElementById('pv-store');
    st.innerHTML = STORE.map(function (s) {
      var n = bytes(s.key); if (s.antigua && n < 0) return '';
      return '<li><div><b>' + esc(s.nombre) + '</b> <code>' + esc(s.key) + '</code><p>' + esc(s.fin) + '</p><small>' + (n < 0 ? 'No guardado en este navegador' : 'Guardado en este navegador · ' + size(n)) + '</small></div>'
        + (n < 0 ? '' : '<button type="button" class="pv-btn pv-mini" data-pv-wipe="' + esc(s.key) + '" aria-label="Borrar ' + esc(s.nombre) + '">Borrar</button>') + '</li>';
    }).join('');
    var on = d.ext.wikimedia, cb = document.getElementById('pv-wiki');
    cb.checked = on;
    var last = d.hist.filter(function (h) { return h.fin === 'wikimedia'; }).pop();
    document.getElementById('pv-wiki-s').textContent = on ? 'Activado' + (last ? ' desde el ' + when(last.t) : '') : 'Desactivado' + (last ? ' desde el ' + when(last.t) : ' (valor predeterminado)');
    document.getElementById('pv-gpc').hidden = !P.gpc();
    var h = document.getElementById('pv-hist');
    h.innerHTML = d.hist.length ? d.hist.slice().reverse().map(function (x) {
      return '<li><b>' + (x.ok ? 'Activaste' : 'Desactivaste') + '</b> ' + esc(PURPOSES[x.fin].nombre.toLowerCase()) + '<small>' + esc(when(x.t)) + ' · ' + esc(METODOS[x.metodo] || x.metodo) + ' · ' + esc(x.doc) + (x.retirado ? ' · retirado el ' + esc(when(x.retirado)) : '') + '</small></li>';
    }).join('') : '<li class="pv-empty">Todavía no has tomado ninguna decisión: todo sigue con los valores predeterminados.</li>';
  },
  say: function (m) { var el = document.getElementById('pv-msg'); if (el) { el.textContent = m; clearTimeout(P._sayT); P._sayT = setTimeout(function () { el.textContent = ''; }, 6000); } },
  open: function () {
    var w = P.build();
    P.from = document.activeElement;
    P.refresh();
    w.hidden = false;
    requestAnimationFrame(function () { w.classList.add('open'); });
    document.documentElement.classList.add('pv-lock');
    setTimeout(function () { var b = w.querySelector('.pv-body'); if (b) b.scrollTop = 0; var c = w.querySelector('.pv-x'); if (c) c.focus({ preventScroll: true }); }, 30);
  },
  close: function () {
    if (!P.el || P.el.hidden) return;
    P.el.classList.remove('open'); P.el.hidden = true;
    document.documentElement.classList.remove('pv-lock');
    var f = P.from;
    if (f && f.closest && f.offsetParent === null) {   // se abrió desde un menú que ya se cerró: el foco vuelve al botón del menú
      var m = f.closest('[role="menu"]');
      f = m ? (document.querySelector('[aria-controls="' + m.id + '"]') || m.parentNode.querySelector('[aria-haspopup]')) : null;
    }
    if (f && f.focus && document.contains(f)) f.focus({ preventScroll: true });
  },
  isOpen: function () { return !!(P.el && !P.el.hidden); },
};
// teclado con el panel abierto: se registra antes que los atajos del simulador (fase de captura) para que Escape cierre
// primero el panel y ningún atajo reaccione mientras se usa; el foco no sale del panel
window.addEventListener('keydown', function (e) {
  if (!P.isOpen()) return;
  var w = P.el;
  if (e.key === 'Escape') { e.preventDefault(); P.close(); }
  else if (e.key === 'Tab') {
    var f = Array.prototype.filter.call(w.querySelectorAll('a[href],button:not([disabled]),input:not([disabled])'), function (x) { return x.offsetParent !== null; });
    if (f.length && (!w.contains(document.activeElement) || (e.shiftKey ? document.activeElement === f[0] : document.activeElement === f[f.length - 1]))) { e.preventDefault(); f[e.shiftKey ? f.length - 1 : 0].focus(); }
  }
  e.stopImmediatePropagation();
}, true);
// cualquier botón o enlace con data-privacy-open abre el panel (en el simulador y en las páginas legales)
document.addEventListener('click', function (e) {
  var o = e.target.closest && e.target.closest('[data-privacy-open]');
  if (!o) return;
  e.preventDefault();
  P.open();
});
// otra pestaña de SOLARIS (por ejemplo, una página legal) cambió las preferencias
window.addEventListener('storage', function (e) { if (e.key === KEY || e.key === null) { emit(); P.refresh(); } });
// enlace directo al panel: …/#privacidad
function hashOpen() { if (location.hash === '#privacidad' || location.hash === '#gestionar-privacidad') P.open(); }
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', hashOpen); else hashOpen();
window.addEventListener('hashchange', hashOpen);
window.SolarisPrivacy = P;
})();
