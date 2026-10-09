(function(){
'use strict';
const DEG = Math.PI / 180, RAD = 180 / Math.PI, TAU = Math.PI * 2;
const AU_KM = 149597870.7;
const J2000 = 2451545.0;
const G_CONST = 6.674e-11;
const clamp = (x, a, b) => (x < a ? a : x > b ? b : x);
const lerp = (a, b, t) => a + (b - a) * t;
const smoothstep = (a, b, x) => { const t = clamp((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); };
const easeInOut = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);
const easeOut = t => 1 - Math.pow(1 - t, 3);
const wrap360 = a => ((a % 360) + 360) % 360;
const wrapPi = a => { a = (a + Math.PI) % TAU; if (a < 0) a += TAU; return a - Math.PI; };
const V = {
add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
scale: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
len: a => Math.hypot(a[0], a[1], a[2]),
norm: a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; },
lerp: (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t],
dist: (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]),
rot: (v, k, ang) => {
const c = Math.cos(ang), s = Math.sin(ang), d = V.dot(k, v), x = V.cross(k, v);
return [v[0] * c + x[0] * s + k[0] * d * (1 - c), v[1] * c + x[1] * s + k[1] * d * (1 - c), v[2] * c + x[2] * s + k[2] * d * (1 - c)];
},
};
const M4 = {
ident() { const m = new Float64Array(16); m[0] = m[5] = m[10] = m[15] = 1; return m; },
mul(a, b, out) {
out = out || new Float64Array(16);
for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
out[c * 4 + r] = a[r] * b[c * 4] + a[4 + r] * b[c * 4 + 1] + a[8 + r] * b[c * 4 + 2] + a[12 + r] * b[c * 4 + 3];
}
return out;
},
perspInf(fovy, aspect, near, out) {
const m = out || new Float64Array(16); m.fill(0);
const f = 1 / Math.tan(fovy / 2);
m[0] = f / aspect; m[5] = f; m[10] = -1; m[11] = -1; m[14] = -2 * near;
return m;
},
viewRot(f, up, out) {
const m = out || new Float64Array(16); m.fill(0);
let r = V.norm(V.cross(f, up));
if (!isFinite(r[0]) || V.len(V.cross(f, up)) < 1e-9) r = V.norm(V.cross(f, [1, 0, 0]));
const u = V.cross(r, f);
m[0] = r[0]; m[4] = r[1]; m[8] = r[2];
m[1] = u[0]; m[5] = u[1]; m[9] = u[2];
m[2] = -f[0]; m[6] = -f[1]; m[10] = -f[2];
m[15] = 1;
return m;
},
fromBasis(x, y, z, s, t, out) {
const m = out || new Float64Array(16);
m[0] = x[0] * s; m[1] = x[1] * s; m[2] = x[2] * s; m[3] = 0;
m[4] = y[0] * s; m[5] = y[1] * s; m[6] = y[2] * s; m[7] = 0;
m[8] = z[0] * s; m[9] = z[1] * s; m[10] = z[2] * s; m[11] = 0;
m[12] = t[0]; m[13] = t[1]; m[14] = t[2]; m[15] = 1;
return m;
},
translate(t, s, out) {
const m = out || new Float64Array(16); m.fill(0);
m[0] = m[5] = m[10] = (s == null ? 1 : s); m[15] = 1; m[12] = t[0]; m[13] = t[1]; m[14] = t[2];
return m;
},
xform(m, p) { // punto
const x = p[0], y = p[1], z = p[2];
const w = m[3] * x + m[7] * y + m[11] * z + m[15];
return [m[0] * x + m[4] * y + m[8] * z + m[12], m[1] * x + m[5] * y + m[9] * z + m[13], m[2] * x + m[6] * y + m[10] * z + m[14], w];
},
};
function makeRng(seed) {
let a = seed >>> 0;
return function () {
a |= 0; a = (a + 0x6D2B79F5) | 0;
let t = Math.imul(a ^ (a >>> 15), 1 | a);
t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
};
}
function gauss(rng) { let u = 0, v = 0; while (u === 0) u = rng(); while (v === 0) v = rng(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * v); }
function hexToRgb(h) { const n = parseInt(h.slice(1), 16); return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]; }
function srgbToLin(c) { return c.map(x => (x <= 0.04045 ? x / 12.92 : Math.pow((x + 0.055) / 1.055, 2.4))); }
const LOCALES = {
es: { name: 'Español', locale: 'es-MX', ui: {}, patterns: [], data: null },
};
const I18N = {
lang: 'es', cache: new Map(), missing: new Set(), busy: false,
L() { return LOCALES[this.lang] || LOCALES.es; },
loc() { return this.L().locale || 'es-MX'; },
t(s, ...a) {
let r = s;
if (this.lang !== 'es') { const u = this.L().ui[s]; if (u != null) r = u; else { const p = this.pattern(s); if (p != null) r = p; } }
return a.length ? r.replace(/\{(\d)\}/g, (m, i) => a[i] != null ? a[i] : '') : r;
},
pattern(s) {
for (const [re, rep] of this.L().patterns || []) { re.lastIndex = 0; if (re.test(s)) { re.lastIndex = 0; return s.replace(re, typeof rep === 'function' ? (...m) => rep(...m) : rep); } }
return null;
},
tr(text) {
if (this.lang === 'es' || !text) return text;
const c = this.cache.get(text); if (c !== undefined) return c;
const m = /^(\s*)([\s\S]*?)(\s*)$/.exec(text), core = m[2];
let out = text;
if (core && /[A-Za-zÁÉÍÓÚáéíóúñÑ]/.test(core)) {
const u = this.L().ui[core];
if (u != null) out = m[1] + u + m[3];
else { const p = this.pattern(core); if (p != null && p !== core) out = m[1] + p + m[3]; else if (this.track) this.missing.add(core); }
}
if (this.cache.size > 6000) this.cache.clear();
this.cache.set(text, out);
return out;
},
ATTRS: ['aria-label', 'placeholder', 'data-tip', 'title', 'alt'],
skip(el) { return el && el.closest && el.closest('[data-no-i18n], script, style, #labels'); },
textNode(n) {
if (this.skip(n.parentElement)) return;
if (n.__o !== undefined && n.data === n.__o) return;          // ya traducido por nosotros
const src = n.data; n.__es = src;
const out = this.tr(src);
n.__o = out; if (out !== src) n.data = out;
},
attr(el, a) {
const v = el.getAttribute(a); if (v == null) return;
el.__ao = el.__ao || {}; el.__ae = el.__ae || {};
if (el.__ao[a] !== undefined && v === el.__ao[a]) return;
el.__ae[a] = v; const out = this.tr(v); el.__ao[a] = out;
if (out !== v) el.setAttribute(a, out);
},
tree(root) {
if (!root) return;
if (root.nodeType === 3) return this.textNode(root);
if (root.nodeType !== 1 || this.skip(root)) return;
const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT); let n;
while ((n = w.nextNode())) this.textNode(n);
const els = [root, ...root.querySelectorAll('*')];
for (const el of els) for (const a of this.ATTRS) if (el.hasAttribute && el.hasAttribute(a)) this.attr(el, a);
},
restore(root) {
const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT); let n;
while ((n = w.nextNode())) if (n.__es !== undefined) { if (n.data !== n.__es && n.data === n.__o) n.data = n.__es; n.__o = undefined; n.__es = undefined; }
root.querySelectorAll('*').forEach(el => { if (el.__ae) { for (const a in el.__ae) if (el.getAttribute(a) === el.__ao[a]) el.setAttribute(a, el.__ae[a]); el.__ae = el.__ao = null; } });
},
observe() {
if (this.obs) return;
this.obs = new MutationObserver(list => {
if (this.lang === 'es' || this.busy) return;
this.busy = true;
for (const m of list) {
if (m.type === 'childList') m.addedNodes.forEach(n => this.tree(n));
else if (m.type === 'characterData') this.textNode(m.target);
else if (m.type === 'attributes' && !this.skip(m.target)) this.attr(m.target, m.attributeName);
}
this.busy = false;
});
this.obs.observe(document.body, { childList: true, subtree: true, characterData: true, attributes: true, attributeFilter: this.ATTRS });
},
swapData(lang) {
const D = (LOCALES[lang] && LOCALES[lang].data) || {};
const keep = (o, keys) => { if (!o.__es) { o.__es = {}; keys.forEach(k => { if (k in o) o.__es[k] = JSON.parse(JSON.stringify(o[k])); }); } };
const put = (o, tr, keys) => { keys.forEach(k => { if (o.__es && k in o.__es) o[k] = JSON.parse(JSON.stringify(o.__es[k])); }); if (tr) for (const k in tr) { if (k === 'info' || k === 'craft') Object.assign(o[k], tr[k]); else o[k] = tr[k]; } };
for (const id in BODY) { const b = BODY[id]; keep(b, ['name', 'short', 'sub', 'info', 'craft']); put(b, D.bodies && D.bodies[id], ['name', 'short', 'sub', 'info', 'craft']);
b.aka = Array.from(new Set([...(b.aka || []), b.__es.name, ...(D.bodies && D.bodies[id] && D.bodies[id].name ? [D.bodies[id].name] : [])])); }
const swapMap = (obj, tr) => { if (!obj) return; if (!obj.__es) Object.defineProperty(obj, '__es', { value: Object.assign({}, obj), enumerable: false }); Object.assign(obj, obj.__es, tr || {}); };
swapMap(TYPE_LABEL, D.typeLabel);
if (typeof DETENT_NAMES !== 'undefined') swapMap(DETENT_NAMES, D.detents);
const swapPairs = (arr, tr) => { if (!arr) return; arr.forEach(p => { if (p.__es === undefined) Object.defineProperty(p, '__es', { value: p[1], enumerable: false }); p[1] = (tr && tr[p[0]]) || p.__es; }); };
if (typeof CAMS !== 'undefined') swapPairs(CAMS, D.cams);
if (typeof STARTS !== 'undefined') swapPairs(STARTS, D.starts);
for (const k in SCALES) { const sc = SCALES[k]; keep(sc, ['label']); put(sc, D.scales && D.scales[k] ? { label: D.scales[k] } : null, ['label']); }
if (typeof SHIPS !== 'undefined') SHIPS.forEach(s => { keep(s, ['name', 'type', 'desc', 'prop']); put(s, D.ships && D.ships[s.id], ['name', 'type', 'desc', 'prop']); });
if (typeof TOURS !== 'undefined') TOURS.forEach(T => { keep(T, ['name', 'title', 'desc']); put(T, D.tours && D.tours[T.id], ['name', 'title', 'desc']);
T.stops.forEach(st => { keep(st, ['group', 'txt']); const tr = {}, X = D.tourText || {}; if (st.__es.group && X[st.__es.group]) tr.group = X[st.__es.group]; if (st.__es.txt && X[st.__es.txt]) tr.txt = X[st.__es.txt]; put(st, tr, ['group', 'txt']); }); });
CONSTELLATIONS.forEach(c => { keep(c, ['name']); put(c, D.const && D.const[c.id] ? { name: D.const[c.id] } : null, ['name']); });
},
set(lang, initial) {
if (!LOCALES[lang]) lang = 'es';
if (lang === this.lang && !initial) return;
if (!TYPE_LABEL.__es) Object.defineProperty(TYPE_LABEL, '__es', { value: Object.assign({}, TYPE_LABEL), enumerable: false });
const prev = this.lang; this.lang = lang; this.cache.clear();
document.documentElement.lang = lang;
this.swapData(lang);
if (lang === 'es') { if (prev !== 'es') this.restore(document.body); }
else { this.observe(); this.busy = true; this.tree(document.body); this.busy = false; }
if (!initial) document.dispatchEvent(new CustomEvent('solaris:lang', { detail: lang }));
},
};
const t = (s, ...a) => I18N.t(s, ...a);
const SRC_PLANET = 'NASA NSSDCA Planetary Fact Sheet · NASA Science';
const SRC_MOON = 'NASA NSSDCA · JPL Solar System Dynamics (satélites)';
const SRC_SMALL = 'JPL Small-Body Database · NASA Science';
const SRC_MOONCOUNT = 'Recuento de lunas: MPC y NASA, marzo de 2026';
const TYPE_LABEL = {
star: 'Estrella', planet: 'Planeta', dwarf: 'Planeta enano', moon: 'Luna', asteroid: 'Asteroide',
comet: 'Cometa', tno: 'Objeto transneptuniano',
};
const BODIES = [
{
id: 'sol', name: 'Sol', aka: ['sun', 'sol'], type: 'star', sub: 'Estrella enana amarilla (G2V)',
R: 695700, mass: 1.989e30, pole: [286.13, 63.87], W0: 84.176, rotH: 609.12, color: '#ffd27a',
vis: { style: 'sun' },
info: {
desc: 'La estrella del Sistema Solar. Concentra el 99.86 % de toda su masa y obtiene su energía de la fusión de hidrógeno en helio en el núcleo.',
age: '≈ 4,600 millones de años', temp: '≈ 5,500 °C en la fotosfera; ≈ 15 millones °C en el núcleo', grav: 274,
day: 'No aplica', rotTxt: '≈ 25 días en el ecuador; ≈ 35 días cerca de los polos',
yearTxt: '≈ 230 millones de años (vuelta alrededor de la Vía Láctea)', vTxt: '≈ 230 km/s alrededor del centro galáctico',
moons: 'No aplica', comp: '≈ 73 % hidrógeno y ≈ 25 % helio (en masa); el resto, oxígeno, carbono, neón, hierro y otros elementos',
atm: 'Cromosfera y corona; la corona supera el millón de grados', tilt: '7.25° respecto a la eclíptica',
feats: ['Ciclo magnético de unos 11 años con manchas, fulguraciones y eyecciones de masa coronal.', 'Rotación diferencial: el ecuador gira más rápido que los polos.', 'El viento solar se extiende más allá de Plutón y forma la heliosfera.'],
facts: ['La luz solar tarda unos 8 minutos y 20 segundos en llegar a la Tierra.', 'Cada segundo convierte unos 4 millones de toneladas de materia en energía.', 'En su volumen cabrían aproximadamente 1.3 millones de Tierras.'],
}, src: 'NASA Sun Fact Sheet · NASA Science',
},
{
id: 'mercurio', name: 'Mercurio', aka: ['mercury'], type: 'planet', parent: 'sol',
R: 2439.7, mass: 3.3011e23, pole: [281.0103, 61.4155], W0: 329.5988, rotH: 1407.6, color: '#b9b3ab',
orbit: { t: 'jpl', el: [0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593], rt: [0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081] },
vis: { style: 'rock', feat: 'mercurio', c1: '#8f8984', c2: '#5d5853', c3: '#b9b2a8', crater: 1.0, seed: 3.1 },
info: {
desc: 'El planeta más pequeño y el más cercano al Sol. Su superficie, cubierta de cráteres, recuerda a la de la Luna.',
age: '≈ 4,500 millones de años', temp: '−173 °C a 427 °C (media ≈ 167 °C)', grav: 3.7,
day: '≈ 176 días terrestres', year: 87.969, v: 47.36, moons: 0, tilt: '0.03°',
comp: 'Rocoso, con un núcleo metálico enorme que ocupa cerca del 85 % de su radio',
atm: 'Exosfera muy tenue de oxígeno, sodio, hidrógeno, helio y potasio',
feats: ['Cuenca Caloris, un impacto de unos 1,550 km de diámetro.', 'Escarpes lobulados que revelan que el planeta se contrajo al enfriarse.', 'Hielo de agua en cráteres polares en sombra permanente.'],
facts: ['Gira tres veces sobre su eje por cada dos vueltas al Sol (resonancia 3:2).', 'Tiene la órbita más excéntrica de los ocho planetas.', 'Lo han visitado Mariner 10 y MESSENGER; BepiColombo (ESA/JAXA) completa su viaje hacia él.'],
}, src: SRC_PLANET,
},
{
id: 'venus', name: 'Venus', aka: ['venus'], type: 'planet', parent: 'sol',
R: 6051.8, mass: 4.8675e24, pole: [272.76, 67.16], W0: 160.20, rotH: -5832.6, color: '#e8cf98',
orbit: { t: 'jpl', el: [0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255], rt: [0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418] },
vis: { style: 'venus', atm: { color: [1.0, 0.86, 0.6], h: 0.06, k: 1.4 } },
info: {
desc: 'Similar a la Tierra en tamaño, pero con un efecto invernadero desbocado que lo convierte en el planeta más caliente.',
age: '≈ 4,500 millones de años', temp: '≈ 464 °C (casi uniforme en todo el planeta)', grav: 8.87,
day: '≈ 116.75 días terrestres', year: 224.701, v: 35.02, moons: 0, tilt: '177.4° (rotación retrógrada)',
comp: 'Rocoso, con una estructura interna parecida a la terrestre',
atm: '≈ 96.5 % CO₂ y 3.5 % N₂; presión en superficie ≈ 92 veces la terrestre; nubes de ácido sulfúrico',
feats: ['Gira en sentido retrógrado: el Sol sale por el oeste.', 'Las nubes altas dan la vuelta al planeta en unos 4 días (superrotación).', 'Superficie volcánica relativamente joven, con miles de volcanes.'],
facts: ['Su día sidéreo (243 días) es más largo que su año (225 días).', 'Es el objeto natural más brillante del cielo nocturno después de la Luna.', 'La sonda soviética Venera 13 sobrevivió 127 minutos en su superficie en 1982.'],
}, src: SRC_PLANET,
},
{
id: 'tierra', name: 'Tierra', aka: ['earth', 'tierra'], type: 'planet', parent: 'sol',
R: 6371.0, mass: 5.9722e24, pole: [0, 90], W0: 190.147, rotH: 23.9344696, color: '#6fa8ff',
orbit: { t: 'jpl', el: [1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0], rt: [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0] },
vis: { style: 'earth', atm: { color: [0.32, 0.58, 1.0], h: 0.035, k: 1.0 }, clouds: true },
info: {
desc: 'Nuestro planeta: el único lugar conocido con vida y con agua líquida estable en la superficie.',
age: '≈ 4,540 millones de años', temp: 'Media ≈ 15 °C (de −89 °C a ≈ 57 °C registrados)', grav: 9.81,
day: '24 h (día solar medio)', year: 365.256, v: 29.78, moons: 1, tilt: '23.44°',
comp: 'Rocosa: núcleo de hierro y níquel, manto de silicatos y corteza',
atm: '≈ 78 % N₂, 21 % O₂, 0.9 % Ar; trazas de CO₂ y vapor de agua',
feats: ['El 71 % de la superficie está cubierta por océanos.', 'Su campo magnético desvía el viento solar y produce auroras.', 'Tectónica de placas activa que recicla la corteza.'],
facts: ['Es el planeta más denso del Sistema Solar (5.51 g/cm³).', 'La inclinación de su eje (23.44°) causa las estaciones.', 'El día se alarga cerca de 2 milisegundos por siglo por las mareas de la Luna.'],
}, src: SRC_PLANET,
},
{
id: 'marte', name: 'Marte', aka: ['mars'], type: 'planet', parent: 'sol',
R: 3389.5, mass: 6.4171e23, pole: [317.681, 52.887], W0: 176.630, rotH: 24.6229, color: '#d9774b',
orbit: { t: 'jpl', el: [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891], rt: [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343] },
vis: { style: 'rock', feat: 'marte', c1: '#b5603a', c2: '#6e3a26', c3: '#d79a6c', crater: 0.35, seed: 7.7, atm: { color: [0.95, 0.62, 0.45], h: 0.02, k: 0.6 } },
info: {
desc: 'El planeta rojo: un mundo frío y desértico con los volcanes y cañones más grandes del Sistema Solar.',
age: '≈ 4,500 millones de años', temp: '−153 °C a 20 °C (media ≈ −65 °C)', grav: 3.71,
day: '24 h 39 min (sol marciano)', year: 686.98, v: 24.07, moons: 2, tilt: '25.19°',
comp: 'Rocoso; superficie rica en óxidos de hierro que le dan su color',
atm: '≈ 95 % CO₂, 2.8 % N₂, 2 % Ar; presión ≈ 0.6 % de la terrestre',
feats: ['Olympus Mons, un volcán de unos 22 km de altura.', 'Valles Marineris, un sistema de cañones de más de 4,000 km.', 'Casquetes polares de hielo de agua y CO₂.'],
facts: ['Tiene estaciones como la Tierra gracias a una inclinación axial parecida.', 'Sus tormentas de polvo pueden cubrir todo el planeta.', 'Rovers como Curiosity y Perseverance exploran su superficie.'],
}, src: SRC_PLANET,
},
{
id: 'jupiter', name: 'Júpiter', aka: ['jupiter'], type: 'planet', parent: 'sol',
R: 69911, mass: 1.8982e27, pole: [268.057, 64.495], W0: 284.95, rotH: 9.925, color: '#d9b48a',
orbit: { t: 'jpl', el: [5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909], rt: [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106] },
vis: { style: 'gas', feat: 'jupiter', atm: { color: [0.85, 0.78, 0.66], h: 0.02, k: 0.5 }, rings: { in: 1.72, out: 1.81, type: 'faint', alpha: 0.045 } },
info: {
desc: 'El gigante gaseoso: más del doble de masivo que todos los demás planetas juntos.',
age: '≈ 4,500 millones de años', temp: '≈ −110 °C en el nivel de 1 bar', grav: 24.79,
day: '≈ 9 h 56 min', year: 4332.59, v: 13.06, moons: 101, tilt: '3.13°',
comp: 'Principalmente hidrógeno y helio; probable núcleo difuso de elementos pesados',
atm: '≈ 90 % H₂ y ≈ 10 % He; nubes de amoníaco, hidrosulfuro de amonio y agua',
feats: ['La Gran Mancha Roja, una tormenta anticiclónica más grande que la Tierra, observada desde hace más de 150 años.', 'Bandas claras (zonas) y oscuras (cinturones) impulsadas por vientos alternos.', 'La magnetosfera más intensa entre los planetas.'],
facts: ['Es el planeta que gira más rápido: su día dura menos de 10 horas.', 'Tiene un sistema de anillos tenue, descubierto por Voyager 1 en 1979.', 'Europa Clipper (NASA) y JUICE (ESA) viajan hacia él para estudiar sus lunas.'],
}, src: SRC_PLANET + ' · ' + SRC_MOONCOUNT,
},
{
id: 'saturno', name: 'Saturno', aka: ['saturn'], type: 'planet', parent: 'sol',
R: 58232, mass: 5.6834e26, pole: [40.589, 83.537], W0: 38.90, rotH: 10.5606, color: '#e6cf98',
orbit: { t: 'jpl', el: [9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448], rt: [-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794] },
vis: { style: 'gas', feat: 'saturno', atm: { color: [0.92, 0.84, 0.62], h: 0.02, k: 0.5 }, rings: { in: 1.11, out: 2.36, type: 'saturn', alpha: 1 } },
info: {
desc: 'El planeta de los anillos: un gigante gaseoso menos denso que el agua.',
age: '≈ 4,500 millones de años', temp: '≈ −140 °C en el nivel de 1 bar', grav: 10.44,
day: '≈ 10 h 34 min', year: 10759.22, v: 9.68, moons: 285, tilt: '26.73°',
comp: 'Hidrógeno y helio, con un núcleo de roca y hielo',
atm: '≈ 96 % H₂ y ≈ 3 % He; trazas de metano y amoníaco',
feats: ['Anillos de hielo de agua de unos 280,000 km de diámetro y, a menudo, sólo decenas de metros de espesor.', 'Un hexágono persistente de nubes en el polo norte.', 'Vientos ecuatoriales de hasta unos 1,800 km/h.'],
facts: ['Su densidad media (0.69 g/cm³) es menor que la del agua.', 'La sonda Cassini lo estudió de 2004 a 2017.', 'Con 285 lunas conocidas en marzo de 2026, es el planeta con más satélites.'],
}, src: SRC_PLANET + ' · ' + SRC_MOONCOUNT,
},
{
id: 'urano', name: 'Urano', aka: ['uranus'], type: 'planet', parent: 'sol',
R: 25362, mass: 8.681e25, pole: [257.311, -15.175], W0: 203.81, rotH: -17.24, color: '#a7dde4',
orbit: { t: 'jpl', el: [19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503], rt: [-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589] },
vis: { style: 'gas', feat: 'urano', atm: { color: [0.62, 0.9, 0.95], h: 0.03, k: 0.8 }, rings: { in: 1.64, out: 2.0, type: 'uranus', alpha: 0.5 } },
info: {
desc: 'Un gigante de hielo que gira "acostado", con su eje inclinado casi 98°.',
age: '≈ 4,500 millones de años', temp: '≈ −195 °C en el nivel de 1 bar (mínimas de −224 °C)', grav: 8.87,
day: '≈ 17 h 14 min', year: 30688.5, v: 6.80, moons: 29, tilt: '97.77°',
comp: 'Mezcla de agua, metano y amoníaco helados sobre un núcleo rocoso',
atm: '≈ 83 % H₂, 15 % He y 2.3 % CH₄',
feats: ['Cada polo recibe unos 42 años seguidos de luz solar y otros 42 de oscuridad.', '13 anillos estrechos y oscuros.', 'El metano absorbe la luz roja y le da su tono verdeazulado.'],
facts: ['Fue el primer planeta descubierto con telescopio (William Herschel, 1781).', 'Sólo lo ha visitado una nave: Voyager 2, en 1986.'],
}, src: SRC_PLANET + ' · ' + SRC_MOONCOUNT,
},
{
id: 'neptuno', name: 'Neptuno', aka: ['neptune'], type: 'planet', parent: 'sol',
R: 24622, mass: 1.02413e26, pole: [299.36, 43.46], W0: 249.978, rotH: 16.11, color: '#5b7ff0',
orbit: { t: 'jpl', el: [30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574], rt: [0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664] },
vis: { style: 'gas', feat: 'neptuno', atm: { color: [0.42, 0.62, 1.0], h: 0.03, k: 0.9 }, rings: { in: 1.66, out: 2.56, type: 'neptune', alpha: 0.4 } },
info: {
desc: 'El planeta más lejano: un gigante de hielo azul con los vientos más rápidos del Sistema Solar.',
age: '≈ 4,500 millones de años', temp: '≈ −200 °C en el nivel de 1 bar', grav: 11.15,
day: '≈ 16 h 7 min', year: 60195, v: 5.43, moons: 16, tilt: '28.32°',
comp: 'Agua, metano y amoníaco helados sobre un núcleo rocoso',
atm: '≈ 80 % H₂, 19 % He y 1.5 % CH₄',
feats: ['Vientos de más de 2,000 km/h.', 'Manchas oscuras: tormentas que aparecen y desaparecen en pocos años.', 'Cinco anillos principales, oscuros y tenues.'],
facts: ['Fue el primer planeta descubierto por predicción matemática (1846).', 'En 2011 completó su primera órbita desde su descubrimiento.', 'Voyager 2 es la única nave que lo ha visitado (1989).'],
}, src: SRC_PLANET + ' · ' + SRC_MOONCOUNT,
},
{
id: 'pluton', name: 'Plutón', aka: ['pluto', 'pluton'], type: 'dwarf', parent: 'sol',
R: 1188.3, mass: 1.303e22, pole: [132.993, -6.163], W0: 302.695, rotH: -153.2928, color: '#d9c2a6',
orbit: { t: 'jpl', el: [39.48211675, 0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684], rt: [-0.00031596, 0.00005170, 0.00004818, 145.20780515, -0.04062942, -0.01183482] },
vis: { style: 'rock', feat: 'pluton', c1: '#cdb79f', c2: '#7a4a33', c3: '#f2ece4', crater: 0.25, seed: 11.3, atm: { color: [0.5, 0.65, 1.0], h: 0.04, k: 0.35 } },
info: {
desc: 'El planeta enano más conocido del cinturón de Kuiper, con glaciares de nitrógeno y montañas de hielo de agua.',
age: '≈ 4,500 millones de años', temp: '≈ −225 °C', grav: 0.62,
day: '≈ 6.4 días terrestres', year: 90560, v: 4.67, moons: 5, tilt: '122.5°',
comp: 'Roca y hielo de agua; superficie de hielos de nitrógeno, metano y monóxido de carbono',
atm: 'Tenue, de nitrógeno, metano y CO; con capas de neblina azulada',
feats: ['Sputnik Planitia, el lóbulo occidental de la gran región en forma de corazón (Tombaugh Regio).', 'Montañas de hielo de agua de varios kilómetros de altura.'],
facts: ['Se consideró planeta desde 1930 hasta 2006.', 'New Horizons lo sobrevoló el 14 de julio de 2015.', 'Durante unos 20 años de cada órbita está más cerca del Sol que Neptuno.'],
}, src: SRC_PLANET + ' · New Horizons (NASA)',
},
{
id: 'ceres', name: 'Ceres', aka: ['ceres'], type: 'dwarf', parent: 'sol',
R: 469.7, mass: 9.383e20, pole: [291.418, 66.764], W0: 170.65, rotH: 9.074, color: '#a7a39c',
orbit: { t: 'kep', a: 2.7675, e: 0.0785, i: 10.587, node: 80.27, peri: 73.6, tp: '2022-12-07' },
vis: { style: 'rock', feat: 'ceres', c1: '#7e7a74', c2: '#57534e', c3: '#9a958d', crater: 0.9, seed: 21.7 },
info: {
desc: 'El objeto más grande del cinturón de asteroides y el único planeta enano del Sistema Solar interior.',
age: '≈ 4,500 millones de años', temp: '≈ −105 °C de media', grav: 0.28,
year: 1680, v: 17.9, moons: 0, tilt: '≈ 4°',
comp: 'Roca y hielo de agua, con sales y arcillas',
atm: 'Exosfera muy tenue y transitoria de vapor de agua',
feats: ['Manchas brillantes de carbonato de sodio en el cráter Occator.', 'Ahuna Mons, un criovolcán de unos 4 km de altura.'],
facts: ['Fue el primer asteroide descubierto (Giuseppe Piazzi, 1801).', 'La sonda Dawn lo orbitó de 2015 a 2018.'],
}, src: SRC_SMALL + ' · Dawn (NASA)',
},
{
id: 'eris', name: 'Eris', aka: ['eris'], type: 'dwarf', parent: 'sol',
R: 1163, mass: 1.6466e22, rotH: 378.9, color: '#e8e4dc',
orbit: { t: 'kep', a: 67.86, e: 0.4361, i: 44.04, node: 35.95, peri: 151.6, fit: { date: '2020-01-01', ra: 25.3, dec: -1.6 } },
vis: { style: 'rock', feat: 'icy', c1: '#e6e2da', c2: '#bdb7ad', c3: '#f6f4f0', crater: 0.3, seed: 31.1 },
info: {
desc: 'Uno de los planetas enanos más masivos, en una órbita muy alargada e inclinada más allá de Neptuno.',
age: '≈ 4,500 millones de años', temp: '≈ −243 °C a −217 °C', year: 204199, moons: 1,
comp: 'Roca y hielos; superficie de metano congelado', atm: 'Probablemente congelada sobre la superficie la mayor parte de su órbita',
feats: ['Uno de los cuerpos más reflectantes del Sistema Solar.', 'Órbita inclinada 44° respecto a la eclíptica.'],
facts: ['Su descubrimiento en 2005 llevó a la UAI a definir "planeta enano" en 2006.', 'Tiene una luna: Disnomia.'],
}, src: SRC_SMALL,
},
{
id: 'haumea', name: 'Haumea', aka: ['haumea'], type: 'dwarf', parent: 'sol',
R: 798, shape: [1.45, 0.65, 1.07], mass: 4.006e21, rotH: 3.9155, color: '#ecebe6',
orbit: { t: 'kep', a: 43.12, e: 0.1958, i: 28.21, node: 121.9, peri: 239.0, fit: { date: '2020-01-01', ra: 213.5, dec: 16.5 } },
vis: { style: 'rock', feat: 'icy', c1: '#ecebe6', c2: '#c9c6bd', c3: '#ffffff', crater: 0.2, seed: 41.5, ringsSmall: true },
info: {
desc: 'Un planeta enano con forma de balón alargado, deformado por su rotación extremadamente rápida.',
age: '≈ 4,500 millones de años', temp: '≈ −241 °C (estimación)', year: 103774, moons: 2,
comp: 'Roca cubierta por hielo de agua cristalino', atm: 'No se ha detectado',
feats: ['Gira en menos de 4 horas: uno de los cuerpos grandes más rápidos.', 'Tiene un anillo, descubierto en 2017.'],
facts: ['Sus lunas se llaman Hiʻiaka y Namaka.', 'Su nombre procede de la diosa hawaiana del parto.'],
}, src: SRC_SMALL,
},
{
id: 'makemake', name: 'Makemake', aka: ['makemake'], type: 'dwarf', parent: 'sol',
R: 715, mass: null, rotH: 22.83, color: '#d9b49a',
orbit: { t: 'kep', a: 45.43, e: 0.1613, i: 28.98, node: 79.62, peri: 294.8, fit: { date: '2020-01-01', ra: 192.5, dec: 27.0 } },
vis: { style: 'rock', feat: 'icy', c1: '#cfa98d', c2: '#9c7360', c3: '#eedfd2', crater: 0.2, seed: 51.9 },
info: {
desc: 'Uno de los objetos más brillantes del cinturón de Kuiper, de superficie rojiza cubierta de hielo de metano.',
age: '≈ 4,500 millones de años', temp: '≈ −239 °C (estimación)', year: 111845, moons: 1,
comp: 'Superficie de hielos de metano y etano', atm: 'No se ha detectado una atmósfera global',
feats: ['Su superficie rojiza se debe a compuestos orgánicos (tolinas).', 'Tiene una pequeña luna oscura, S/2015 (136472) 1.'],
facts: ['Se descubrió poco después de la Pascua de 2005; su nombre procede de un dios creador de Rapa Nui.'],
}, src: SRC_SMALL,
},
{
id: 'sedna', name: 'Sedna', aka: ['sedna'], type: 'tno', parent: 'sol',
R: 498, mass: null, rotH: 10.27, color: '#d0745a',
orbit: { t: 'kep', a: 506, e: 0.8496, i: 11.93, node: 144.4, peri: 311.3, fit: { date: '2020-01-01', ra: 58.0, dec: 7.6 } },
vis: { style: 'rock', feat: 'icy', c1: '#a8503a', c2: '#6f3022', c3: '#c97a5e', crater: 0.3, seed: 61.2 },
info: {
desc: 'Un objeto lejano y muy rojizo cuya órbita lo lleva a cientos de unidades astronómicas del Sol.',
age: '≈ 4,500 millones de años', temp: '≈ −240 °C (estimación)', yearTxt: '≈ 11,400 años', moons: 0,
comp: 'Roca y hielos; superficie muy rojiza por compuestos orgánicos (tolinas)', atm: 'No disponible',
feats: ['Perihelio de unas 76 UA y afelio de unas 900 UA.', 'Uno de los objetos más rojos del Sistema Solar.'],
facts: ['Llegará a su perihelio hacia 2076.', 'Su órbita podría deberse al paso de una estrella cercana en el pasado o a un planeta lejano aún no descubierto (hipótesis).'],
}, src: SRC_SMALL,
},
{
id: 'quaoar', name: 'Quaoar', aka: ['quaoar'], type: 'tno', parent: 'sol',
R: 555, mass: 1.2e21, rotH: 17.68, color: '#b88b72',
orbit: { t: 'kep', a: 43.69, e: 0.0414, i: 7.99, node: 188.9, peri: 147.5, fit: { date: '2020-01-01', ra: 261.0, dec: -15.5 } },
vis: { style: 'rock', feat: 'icy', c1: '#9c705c', c2: '#6a4a3c', c3: '#c09a86', crater: 0.4, seed: 71.4 },
info: {
desc: 'Un gran objeto del cinturón de Kuiper con anillos inesperadamente lejanos.',
age: '≈ 4,500 millones de años', temp: '≈ −230 °C (estimación)', year: 105490, moons: 1,
comp: 'Roca y hielo de agua cristalino', atm: 'No disponible',
feats: ['Sus anillos están mucho más allá de su límite de Roche, algo inesperado.', 'Tiene una luna llamada Weywot.'],
facts: ['Su nombre procede de una deidad creadora del pueblo tongva.'],
}, src: SRC_SMALL,
},
{
id: 'arrokoth', name: 'Arrokoth', aka: ['arrokoth', 'ultima thule', '2014 mu69'], type: 'tno', parent: 'sol',
R: 9, shape: [1.9, 0.55, 1.05], mass: null, rotH: 15.92, color: '#b56a4f',
orbit: { t: 'kep', a: 44.58, e: 0.0417, i: 2.45, node: 158.9, peri: 174.0, fit: { date: '2019-01-01', ra: 290.0, dec: -20.6 } },
vis: { style: 'rock', feat: 'bilobe', c1: '#a45a40', c2: '#6a3426', c3: '#c4876a', crater: 0.3, seed: 81.0, irregular: 0.12 },
info: {
desc: 'Un pequeño binario de contacto del cinturón de Kuiper: el objeto más lejano visitado por una nave.',
age: '≈ 4,500 millones de años', temp: '≈ −233 °C (estimación)', year: 108843, moons: 0,
comp: 'Hielos de metano y compuestos orgánicos rojizos', atm: 'Ninguna',
feats: ['Dos lóbulos unidos suavemente, de unos 36 km de largo en total.', 'Superficie poco craterizada.'],
facts: ['New Horizons lo sobrevoló el 1 de enero de 2019.', 'Es un planetesimal casi intacto desde la formación del Sistema Solar.'],
}, src: SRC_SMALL + ' · New Horizons (NASA)',
},
{
id: 'vesta', name: 'Vesta', aka: ['vesta', '4 vesta'], type: 'asteroid', parent: 'sol',
R: 262.7, shape: [1.09, 0.85, 1.05], mass: 2.59e20, pole: [309.031, 42.235], W0: 285.39, rotH: 5.342, color: '#bcb6a8',
orbit: { t: 'kep', a: 2.3615, e: 0.0887, i: 7.142, node: 103.80, peri: 151.2, tp: '2021-12-26' },
vis: { style: 'rock', feat: 'vesta', c1: '#9c958a', c2: '#6c665d', c3: '#c4beb1', crater: 0.9, seed: 91.3, irregular: 0.08 },
info: {
desc: 'El segundo objeto más masivo del cinturón de asteroides, un protoplaneta diferenciado.',
age: '≈ 4,560 millones de años', temp: '≈ −190 °C a −20 °C (aprox.)', grav: 0.25, year: 1325.8, v: 19.34, moons: 0,
comp: 'Diferenciado: núcleo de hierro, manto y corteza basáltica', atm: 'Ninguna',
feats: ['Cráter Rheasilvia, de unos 500 km, con un pico central de unos 22 km.', 'Origen de los meteoritos HED que caen a la Tierra.'],
facts: ['Es el asteroide más brillante visto desde la Tierra.', 'La sonda Dawn lo orbitó en 2011 y 2012.'],
}, src: SRC_SMALL + ' · Dawn (NASA)',
},
{
id: 'palas', name: 'Palas', aka: ['pallas', '2 pallas'], type: 'asteroid', parent: 'sol',
R: 256, shape: [1.08, 0.9, 1.02], mass: 2.04e20, rotH: 7.813, color: '#a39e95',
orbit: { t: 'kep', a: 2.7733, e: 0.2302, i: 34.84, node: 172.92, peri: 310.9, tp: '2022-08-24' },
vis: { style: 'rock', feat: 'generic', c1: '#7f7a73', c2: '#55514c', c3: '#9b968e', crater: 0.9, seed: 101.7, irregular: 0.07 },
info: {
desc: 'Uno de los mayores asteroides, con una órbita muy inclinada respecto a la eclíptica.',
age: '≈ 4,500 millones de años', temp: 'No disponible', year: 1686, v: 17.65, moons: 0,
comp: 'Silicatos ricos en carbono', atm: 'Ninguna',
feats: ['Órbita inclinada unos 35°, lo que dificulta visitarlo.', 'Superficie muy craterizada (observaciones del VLT).'],
facts: ['Descubierto por Heinrich Olbers en 1802.'],
}, src: SRC_SMALL,
},
{
id: 'higia', name: 'Higía', aka: ['hygiea', 'higia', '10 hygiea'], type: 'asteroid', parent: 'sol',
R: 217, mass: 8.7e19, rotH: 13.83, color: '#8e8a84',
orbit: { t: 'kep', a: 3.1415, e: 0.1125, i: 3.832, node: 283.2, peri: 312.3, tp: '2023-05-01' },
vis: { style: 'rock', feat: 'generic', c1: '#5f5c58', c2: '#43413e', c3: '#77736e', crater: 0.7, seed: 111.1 },
info: {
desc: 'El cuarto objeto más grande del cinturón de asteroides, casi esférico.',
age: '≈ 4,500 millones de años', temp: 'No disponible', year: 2033, moons: 0,
comp: 'Carbonáceo (tipo C)', atm: 'Ninguna',
feats: ['Su forma casi esférica lo convierte en candidato a planeta enano.', 'Es el mayor miembro de una familia de asteroides.'],
facts: ['Descubierto por Annibale de Gasparis en 1849.'],
}, src: SRC_SMALL,
},
{
id: 'psique', name: 'Psique', aka: ['psyche', '16 psyche'], type: 'asteroid', parent: 'sol',
R: 111, shape: [1.25, 0.85, 1.0], mass: 2.29e19, rotH: 4.196, color: '#b0aaa0',
orbit: { t: 'kep', a: 2.924, e: 0.134, i: 3.095, node: 150.03, peri: 229.3, tp: '2021-12-01' },
vis: { style: 'rock', feat: 'metal', c1: '#8d877e', c2: '#5e5951', c3: '#c2bcb0', crater: 0.7, seed: 121.5, irregular: 0.12 },
info: {
desc: 'Un asteroide que podría ser el núcleo metálico expuesto de un protoplaneta.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Probablemente rico en metales (hierro y níquel) mezclados con roca', atm: 'Ninguna',
feats: ['Uno de los objetos metálicos más grandes conocidos.', 'La sonda Psyche de la NASA tiene prevista su llegada en 2029.'],
facts: ['Descubierto por Annibale de Gasparis en 1852.'],
}, src: SRC_SMALL,
},
{
id: 'eros', name: 'Eros', aka: ['eros', '433 eros'], type: 'asteroid', parent: 'sol',
R: 8.42, shape: [2.0, 0.66, 0.66], mass: 6.687e15, rotH: 5.270, color: '#b59a7d',
orbit: { t: 'kep', a: 1.4580, e: 0.2227, i: 10.83, node: 304.30, peri: 178.88, fit: { date: '2012-01-31', earth: true } },
vis: { style: 'rock', feat: 'generic', c1: '#9c8269', c2: '#6c5845', c3: '#bba184', crater: 0.6, seed: 131.9, irregular: 0.22 },
info: {
desc: 'Asteroide cercano a la Tierra con forma alargada; el primero en el que se posó una nave.',
age: '≈ 4,500 millones de años', temp: 'No disponible', year: 643.2, v: 24.36, moons: 0,
comp: 'Silicatos (tipo S)', atm: 'Ninguna',
feats: ['Mide unos 34 × 11 × 11 km.', 'Superficie cubierta de regolito y rocas sueltas.'],
facts: ['NEAR Shoemaker aterrizó en él en febrero de 2001.', 'Su posición se ajusta aquí a su acercamiento a la Tierra de enero de 2012 (aproximación).'],
}, src: SRC_SMALL,
},
{
id: 'bennu', name: 'Bennu', aka: ['bennu', '101955 bennu'], type: 'asteroid', parent: 'sol',
R: 0.245, shape: [1.05, 0.95, 1.05], mass: 7.329e10, rotH: 4.296, color: '#6f6b67',
orbit: { t: 'kep', a: 1.1264, e: 0.2037, i: 6.035, node: 2.061, peri: 66.22, fit: { date: '2005-09-20', earth: true } },
vis: { style: 'rock', feat: 'rubble', c1: '#4b4845', c2: '#2f2d2b', c3: '#6a6662', crater: 0.4, seed: 141.2, irregular: 0.1 },
info: {
desc: 'Un pequeño asteroide carbonáceo cercano a la Tierra, del que se trajeron muestras.',
age: '≈ 4,500 millones de años (material)', temp: 'No disponible', year: 436.6, moons: 0,
comp: 'Carbonáceo (tipo B), con minerales hidratados y compuestos orgánicos', atm: 'Ninguna',
feats: ['Forma de trompo y superficie cubierta de rocas.', 'Expulsa pequeñas partículas al espacio.'],
facts: ['OSIRIS-REx trajo 121.6 g de muestra a la Tierra en septiembre de 2023.', 'Tiene una probabilidad pequeña de impactar la Tierra a finales del siglo XXII.'],
}, src: SRC_SMALL + ' · OSIRIS-REx (NASA)',
},
{
id: 'ryugu', name: 'Ryugu', aka: ['ryugu', '162173 ryugu'], type: 'asteroid', parent: 'sol',
R: 0.45, shape: [1.03, 0.94, 1.03], mass: 4.5e11, rotH: 7.627, color: '#6a6561',
orbit: { t: 'kep', a: 1.1896, e: 0.1902, i: 5.884, node: 251.59, peri: 211.43, tp: '2021-11-26' },
vis: { style: 'rock', feat: 'rubble', c1: '#47433f', c2: '#2c2a28', c3: '#625d58', crater: 0.4, seed: 151.6, irregular: 0.1 },
info: {
desc: 'Asteroide carbonáceo con forma de trompo visitado por la misión japonesa Hayabusa2.',
age: '≈ 4,500 millones de años (material)', temp: 'No disponible', year: 474, moons: 0,
comp: 'Carbonáceo (tipo C)', atm: 'Ninguna',
feats: ['Forma de trompo con una cresta ecuatorial.', 'Muy poroso: es un "montón de escombros".'],
facts: ['Hayabusa2 (JAXA) trajo sus muestras a la Tierra en diciembre de 2020.'],
}, src: SRC_SMALL + ' · Hayabusa2 (JAXA)',
},
{
id: 'apofis', name: 'Apofis', aka: ['apophis', '99942 apophis'], type: 'asteroid', parent: 'sol',
R: 0.17, shape: [1.25, 0.85, 1.0], mass: null, rotH: 30.56, color: '#9e8f80',
orbit: { t: 'kep', a: 0.9224, e: 0.1915, i: 3.339, node: 203.96, peri: 126.60, fit: { date: '2029-04-13', earth: true } },
vis: { style: 'rock', feat: 'generic', c1: '#8a7b6c', c2: '#5d5249', c3: '#a8998a', crater: 0.4, seed: 161.3, irregular: 0.15 },
info: {
desc: 'Asteroide cercano a la Tierra que en 2029 pasará excepcionalmente cerca de nuestro planeta.',
age: '≈ 4,500 millones de años (material)', temp: 'No disponible', year: 323.6, moons: 0,
comp: 'Silicatos (tipo Sq)', atm: 'Ninguna',
feats: ['El 13 de abril de 2029 pasará a unos 32,000 km de la superficie terrestre, más cerca que los satélites geoestacionarios.', 'Gira de forma caótica ("tumbling"), sin un eje fijo.'],
facts: ['Los cálculos actuales descartan un impacto durante al menos los próximos 100 años.', 'OSIRIS-APEX lo estudiará tras su paso de 2029.'],
}, src: SRC_SMALL,
},
{
id: 'halley', name: '1P/Halley', short: 'Halley', aka: ['halley', 'cometa halley'], type: 'comet', parent: 'sol',
R: 5.5, shape: [1.4, 0.75, 0.75], mass: 2.2e14, rotH: 52.8, color: '#bfe3ff',
orbit: { t: 'kep', a: 17.834, e: 0.96714, i: 162.262, node: 58.420, peri: 111.332, tp: '1986-02-09' },
vis: { style: 'rock', feat: 'comet', c1: '#3a3632', c2: '#22201e', c3: '#57524c', crater: 0.2, seed: 171.5, irregular: 0.22 },
info: {
desc: 'El cometa periódico más famoso: regresa al Sistema Solar interior cada 75 o 76 años.',
age: '≈ 4,600 millones de años (material primitivo)', temp: 'No disponible', rotTxt: '≈ 2.2 días (rotación compleja)', yearTxt: '≈ 75–76 años', moons: 0,
comp: 'Hielo de agua, CO, CO₂ y polvo oscuro', atm: 'Coma y colas de gas y polvo cuando se acerca al Sol',
feats: ['Núcleo muy oscuro: refleja sólo ≈ 4 % de la luz.', 'Órbita retrógrada inclinada 162°.'],
facts: ['Último perihelio: 9 de febrero de 1986; el próximo será hacia julio de 2061.', 'Origina las lluvias de meteoros Eta Acuáridas y Oriónidas.', 'La sonda Giotto fotografió su núcleo en 1986.'],
}, src: SRC_SMALL,
},
{
id: 'halebopp', name: 'C/1995 O1 Hale-Bopp', short: 'Hale-Bopp', aka: ['hale-bopp', 'hale bopp'], type: 'comet', parent: 'sol',
R: 30, shape: [1.1, 0.9, 1.0], mass: null, rotH: 11.35, color: '#d8f0ff',
orbit: { t: 'kep', a: 186.5, e: 0.99510, i: 89.43, node: 282.47, peri: 130.59, tp: '1997-04-01' },
vis: { style: 'rock', feat: 'comet', c1: '#3a3632', c2: '#22201e', c3: '#57524c', crater: 0.2, seed: 181.1, irregular: 0.12 },
info: {
desc: 'Uno de los cometas más brillantes del siglo XX, visible a simple vista durante unos 18 meses.',
age: '≈ 4,600 millones de años (material primitivo)', temp: 'No disponible', yearTxt: '≈ 2,400 años', moons: 0,
comp: 'Hielos y polvo; núcleo de unos 60 km (estimación)', atm: 'Coma y colas de gas, polvo y sodio',
feats: ['Uno de los núcleos cometarios más grandes conocidos.', 'Mostró una tercera cola formada por sodio.'],
facts: ['Perihelio: 1 de abril de 1997, a 0.91 UA del Sol.'],
}, src: SRC_SMALL,
},
{
id: 'c67p', name: '67P/Churyumov-Gerasimenko', short: '67P', aka: ['67p', 'churyumov', 'gerasimenko', 'rosetta'], type: 'comet', parent: 'sol',
R: 2.0, shape: [1.3, 0.85, 0.8], mass: 1.0e13, rotH: 12.4, color: '#cfe8ff',
orbit: { t: 'kep', a: 3.4630, e: 0.6410, i: 7.040, node: 50.136, peri: 12.780, tp: '2021-11-02' },
vis: { style: 'rock', feat: 'bilobe', c1: '#3c3936', c2: '#252321', c3: '#59554f', crater: 0.15, seed: 191.7, irregular: 0.1 },
info: {
desc: 'Cometa de periodo corto explorado de cerca por la misión Rosetta de la ESA.',
age: '≈ 4,600 millones de años (material primitivo)', temp: 'No disponible', yearTxt: '≈ 6.44 años', moons: 0,
comp: 'Hielo y polvo muy poroso, con abundantes compuestos orgánicos', atm: 'Coma de gas y polvo cerca del perihelio',
feats: ['Forma de "patito de goma" con dos lóbulos.', 'Chorros de gas y polvo al acercarse al Sol.'],
facts: ['Rosetta lo orbitó de 2014 a 2016; el módulo Philae aterrizó el 12 de noviembre de 2014.'],
}, src: SRC_SMALL + ' · Rosetta (ESA)',
},
{
id: 'encke', name: '2P/Encke', short: 'Encke', aka: ['encke'], type: 'comet', parent: 'sol',
R: 2.4, mass: null, rotH: 11.1, color: '#cfe8ff',
orbit: { t: 'kep', a: 2.2178, e: 0.8471, i: 11.78, node: 334.57, peri: 186.55, tp: '2023-10-22' },
vis: { style: 'rock', feat: 'comet', c1: '#3a3632', c2: '#22201e', c3: '#57524c', crater: 0.2, seed: 201.2, irregular: 0.15 },
info: {
desc: 'El cometa con el periodo orbital más corto entre los cometas brillantes conocidos.',
age: '≈ 4,600 millones de años (material primitivo)', temp: 'No disponible', yearTxt: '≈ 3.3 años', moons: 0,
comp: 'Hielos y polvo', atm: 'Coma cerca del perihelio',
feats: ['Origina las lluvias de meteoros Táuridas.', 'Su núcleo mide unos 4.8 km.'],
facts: ['Fue el segundo cometa cuyo regreso se predijo, después del Halley.'],
}, src: SRC_SMALL,
},
{
id: 'swifttuttle', name: '109P/Swift-Tuttle', short: 'Swift-Tuttle', aka: ['swift-tuttle', 'perseidas'], type: 'comet', parent: 'sol',
R: 13, mass: null, rotH: null, color: '#cfe8ff',
orbit: { t: 'kep', a: 26.092, e: 0.9632, i: 113.45, node: 139.38, peri: 152.98, tp: '1992-12-11' },
vis: { style: 'rock', feat: 'comet', c1: '#3a3632', c2: '#22201e', c3: '#57524c', crater: 0.2, seed: 211.8, irregular: 0.12 },
info: {
desc: 'Cometa progenitor de las Perseidas, la lluvia de estrellas de agosto.',
age: '≈ 4,600 millones de años (material primitivo)', temp: 'No disponible', yearTxt: '≈ 133 años', moons: 0,
comp: 'Hielos y polvo; núcleo de unos 26 km', atm: 'Coma y cola cerca del perihelio',
feats: ['El objeto conocido más grande que pasa repetidamente cerca de la órbita terrestre.', 'Origina las Perseidas cada mes de agosto.'],
facts: ['Último perihelio en diciembre de 1992; el próximo será en 2126.'],
}, src: SRC_SMALL,
},
{
id: 'luna', name: 'Luna', aka: ['moon', 'luna'], type: 'moon', parent: 'tierra',
R: 1737.4, mass: 7.342e22, locked: true, color: '#cfcac2', orbit: { t: 'luna', a: 384400, P: 27.321661 },
vis: { style: 'rock', feat: 'luna', c1: '#a7a29a', c2: '#57534e', c3: '#d6d1c8', crater: 1.0, seed: 1.7 },
info: {
desc: 'El único satélite natural de la Tierra y el único cuerpo fuera de ella pisado por seres humanos.',
age: '≈ 4,500 millones de años', temp: '−173 °C a 127 °C', grav: 1.62, day: '≈ 29.5 días (mes sinódico)', v: 1.022, moons: 0,
comp: 'Rocosa, con un núcleo metálico pequeño', atm: 'Exosfera extremadamente tenue', tilt: '6.7° respecto a su órbita',
feats: ['Rotación síncrona: siempre muestra la misma cara a la Tierra.', 'Los mares son llanuras de basalto formadas por antiguas erupciones.', 'Cráteres jóvenes como Tycho tienen rayos brillantes de material eyectado.'],
facts: ['Se aleja de la Tierra unos 3.8 cm por año.', 'Doce astronautas caminaron sobre ella entre 1969 y 1972.', 'Probablemente se formó tras el choque de un cuerpo del tamaño de Marte con la Tierra primitiva.'],
}, src: SRC_MOON,
},
{
id: 'fobos', name: 'Fobos', aka: ['phobos'], type: 'moon', parent: 'marte',
R: 11.08, shape: [1.22, 0.82, 1.0], mass: 1.0659e16, locked: true, color: '#8d8378',
orbit: { t: 'moon', a: 9376, P: 0.31891023, i: 1.08, node: 0, L0: 35 },
vis: { style: 'rock', feat: 'generic', c1: '#6e665e', c2: '#4a443f', c3: '#8a8178', crater: 0.8, seed: 221.4, irregular: 0.12 },
info: {
desc: 'La mayor y más cercana luna de Marte: un cuerpo irregular cubierto de polvo y cráteres.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Probablemente roca rica en carbono, muy porosa', atm: 'Ninguna',
feats: ['Cráter Stickney, de unos 9 km de diámetro.', 'Estrías paralelas que recorren su superficie.'],
facts: ['Orbita más rápido de lo que Marte gira: desde Marte sale por el oeste dos veces al día.', 'Se acerca a Marte cerca de 1.8 m por siglo; podría fragmentarse en un anillo en unos 50 millones de años.'],
}, src: SRC_MOON,
},
{
id: 'deimos', name: 'Deimos', aka: ['deimos'], type: 'moon', parent: 'marte',
R: 6.2, shape: [1.21, 0.89, 0.92], mass: 1.4762e15, locked: true, color: '#9a8f82',
orbit: { t: 'moon', a: 23463, P: 1.263, i: 1.79, node: 0, L0: 200 },
vis: { style: 'rock', feat: 'generic', c1: '#7b7268', c2: '#58514a', c3: '#978d82', crater: 0.4, seed: 231.3, irregular: 0.1 },
info: {
desc: 'La pequeña luna exterior de Marte, de superficie suavizada por el polvo.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Probablemente roca rica en carbono', atm: 'Ninguna',
feats: ['Superficie cubierta por una capa de regolito de decenas de metros.', 'Mide unos 15 × 12 × 11 km.'],
facts: ['Desde Marte se ve como una estrella brillante.', 'Su nombre significa "terror" en griego.'],
}, src: SRC_MOON,
},
{
id: 'io', name: 'Ío', aka: ['io'], type: 'moon', parent: 'jupiter',
R: 1821.6, mass: 8.9319e22, locked: true, color: '#e8d27a',
orbit: { t: 'moon', a: 421700, P: 1.769137786, i: 0.04, node: 0, L0: 106.07719 },
vis: { style: 'rock', feat: 'io', c1: '#d3bd70', c2: '#8a5a2a', c3: '#efe6cc', crater: 0.0, seed: 241.6 },
info: {
desc: 'El cuerpo con mayor actividad volcánica del Sistema Solar.',
age: '≈ 4,500 millones de años', temp: '≈ −150 °C de media; algunas lavas superan 1,200 °C', grav: 1.796, moons: 0,
comp: 'Silicatos con núcleo de hierro; superficie de azufre y compuestos sulfurosos', atm: 'Muy tenue, de dióxido de azufre',
feats: ['Más de 400 volcanes activos.', 'Las mareas de Júpiter calientan su interior.'],
facts: ['Algunas columnas volcánicas alcanzan cientos de kilómetros de altura.', 'Galileo Galilei la descubrió en 1610 junto con Europa, Ganímedes y Calisto.'],
}, src: SRC_MOON,
},
{
id: 'europa', name: 'Europa', aka: ['europa'], type: 'moon', parent: 'jupiter',
R: 1560.8, mass: 4.7998e22, locked: true, color: '#e9dfcf',
orbit: { t: 'moon', a: 671034, P: 3.551181, i: 0.47, node: 0, L0: 175.73161 },
vis: { style: 'rock', feat: 'europa', c1: '#e7dccb', c2: '#9a5a3a', c3: '#f5efe6', crater: 0.05, seed: 251.2 },
info: {
desc: 'Una luna helada con un probable océano global bajo su corteza: uno de los lugares más prometedores para buscar vida.',
age: '≈ 4,500 millones de años', temp: '≈ −160 °C en el ecuador; −220 °C en los polos', grav: 1.315, moons: 0,
comp: 'Corteza de hielo de agua sobre un océano salado y un interior rocoso', atm: 'Exosfera tenue de oxígeno',
feats: ['Superficie muy lisa cruzada por grietas y crestas rojizas.', 'Océano bajo una capa de hielo de entre 15 y 25 km (estimación).'],
facts: ['Su océano podría contener el doble de agua que todos los océanos terrestres.', 'Es el objetivo principal de la misión Europa Clipper.'],
}, src: SRC_MOON,
},
{
id: 'ganimedes', name: 'Ganímedes', aka: ['ganymede', 'ganimedes'], type: 'moon', parent: 'jupiter',
R: 2634.1, mass: 1.4819e23, locked: true, color: '#b3a796',
orbit: { t: 'moon', a: 1070412, P: 7.15455296, i: 0.18, node: 0, L0: 120.55883 },
vis: { style: 'rock', feat: 'ganimedes', c1: '#8f8475', c2: '#5d5448', c3: '#cfc6b8', crater: 0.6, seed: 261.9 },
info: {
desc: 'La luna más grande del Sistema Solar, mayor incluso que Mercurio.',
age: '≈ 4,500 millones de años', temp: '≈ −160 °C (aprox.)', grav: 1.428, moons: 0,
comp: 'Roca y hielo; núcleo de hierro líquido', atm: 'Exosfera tenue de oxígeno',
feats: ['Regiones oscuras antiguas y terrenos claros surcados.', 'Probable océano salado bajo la corteza.'],
facts: ['Es la única luna conocida con campo magnético propio.', 'La misión JUICE (ESA) la orbitará en la década de 2030.'],
}, src: SRC_MOON,
},
{
id: 'calisto', name: 'Calisto', aka: ['callisto', 'calisto'], type: 'moon', parent: 'jupiter',
R: 2410.3, mass: 1.0759e23, locked: true, color: '#7e7468',
orbit: { t: 'moon', a: 1882709, P: 16.6890184, i: 0.19, node: 0, L0: 84.44459 },
vis: { style: 'rock', feat: 'calisto', c1: '#5a5148', c2: '#3a342e', c3: '#c9c0b2', crater: 1.0, seed: 271.4 },
info: {
desc: 'Una de las superficies más craterizadas del Sistema Solar, casi sin cambios en miles de millones de años.',
age: '≈ 4,500 millones de años', temp: '≈ −139 °C', grav: 1.235, moons: 0,
comp: 'Roca y hielo poco diferenciados', atm: 'Exosfera muy tenue de CO₂',
feats: ['Valhalla, una cuenca de impacto con anillos de unos 3,800 km.', 'Posible océano profundo bajo la corteza.'],
facts: ['Orbita fuera de la zona más intensa de radiación de Júpiter.'],
}, src: SRC_MOON,
},
{
id: 'mimas', name: 'Mimas', aka: ['mimas'], type: 'moon', parent: 'saturno',
R: 198.2, mass: 3.7493e19, locked: true, color: '#c8c4bd',
orbit: { t: 'moon', a: 185539, P: 0.942422, i: 1.574, node: 0, L0: 14 },
vis: { style: 'rock', feat: 'mimas', c1: '#bdb8b0', c2: '#8c8780', c3: '#d8d4cc', crater: 1.0, seed: 281.7 },
info: {
desc: 'Una pequeña luna helada dominada por un cráter gigantesco.',
age: '≈ 4,500 millones de años', temp: '≈ −200 °C (aprox.)', moons: 0,
comp: 'Principalmente hielo de agua', atm: 'Ninguna',
feats: ['Cráter Herschel, de 130 km: casi un tercio de su diámetro.', 'Indicios de un océano interno joven (datos de Cassini, 2024).'],
facts: ['Su aspecto recuerda a la "Estrella de la Muerte".', 'Su gravedad abre la división de Cassini en los anillos (resonancia 2:1).'],
}, src: SRC_MOON,
},
{
id: 'encelado', name: 'Encélado', aka: ['enceladus', 'encelado'], type: 'moon', parent: 'saturno',
R: 252.1, mass: 1.0802e20, locked: true, color: '#f2f5f7',
orbit: { t: 'moon', a: 237948, P: 1.370218, i: 0.009, node: 0, L0: 210 },
vis: { style: 'rock', feat: 'encelado', c1: '#eef2f4', c2: '#b9c4cc', c3: '#ffffff', crater: 0.3, seed: 291.2 },
info: {
desc: 'Una luna brillante que expulsa chorros de agua desde un océano oculto.',
age: '≈ 4,500 millones de años', temp: '≈ −201 °C', moons: 0,
comp: 'Hielo de agua sobre un océano global y un núcleo rocoso', atm: 'Penachos de vapor de agua en el polo sur',
feats: ['Géiseres que surgen de fracturas llamadas "rayas de tigre".', 'Su material alimenta el anillo E de Saturno.'],
facts: ['Refleja casi toda la luz que recibe: es uno de los cuerpos más brillantes del Sistema Solar.', 'Cassini detectó en sus penachos sales, compuestos orgánicos y fósforo.'],
}, src: SRC_MOON,
},
{
id: 'tetis', name: 'Tetis', aka: ['tethys', 'tetis'], type: 'moon', parent: 'saturno',
R: 531.1, mass: 6.1745e20, locked: true, color: '#dcdad6',
orbit: { t: 'moon', a: 294619, P: 1.887802, i: 1.12, node: 0, L0: 320 },
vis: { style: 'rock', feat: 'icy', c1: '#d6d3cd', c2: '#a9a59e', c3: '#ecebe7', crater: 0.9, seed: 301.6 },
info: {
desc: 'Una luna de hielo casi puro con un enorme valle y un gran cráter.',
age: '≈ 4,500 millones de años', temp: '≈ −187 °C (aprox.)', moons: 0,
comp: 'Casi completamente hielo de agua', atm: 'Ninguna',
feats: ['Ithaca Chasma, un valle de unos 2,000 km de longitud.', 'Cráter Odysseus, de unos 445 km.'],
facts: ['Comparte órbita con dos pequeñas lunas troyanas: Telesto y Calipso.'],
}, src: SRC_MOON,
},
{
id: 'dione', name: 'Dione', aka: ['dione'], type: 'moon', parent: 'saturno',
R: 561.4, mass: 1.0955e21, locked: true, color: '#d2cfca',
orbit: { t: 'moon', a: 377396, P: 2.736915, i: 0.019, node: 0, L0: 95 },
vis: { style: 'rock', feat: 'icy', c1: '#c9c5bf', c2: '#8f8b85', c3: '#e8e6e2', crater: 0.8, seed: 311.3 },
info: {
desc: 'Una luna de hielo y roca con acantilados brillantes en su hemisferio posterior.',
age: '≈ 4,500 millones de años', temp: '≈ −186 °C (aprox.)', moons: 0,
comp: 'Hielo con una proporción importante de roca', atm: 'Exosfera muy tenue de oxígeno',
feats: ['Acantilados de hielo de cientos de metros de altura.', 'Posible océano interno profundo.'],
facts: ['Comparte órbita con las lunas troyanas Helena y Pólux.'],
}, src: SRC_MOON,
},
{
id: 'rea', name: 'Rea', aka: ['rhea', 'rea'], type: 'moon', parent: 'saturno',
R: 763.8, mass: 2.3065e21, locked: true, color: '#cfccc6',
orbit: { t: 'moon', a: 527108, P: 4.518212, i: 0.345, node: 0, L0: 250 },
vis: { style: 'rock', feat: 'icy', c1: '#c4c0ba', c2: '#8a8680', c3: '#e2e0dc', crater: 1.0, seed: 321.9 },
info: {
desc: 'La segunda luna más grande de Saturno, helada y muy craterizada.',
age: '≈ 4,500 millones de años', temp: '−174 °C a −220 °C', moons: 0,
comp: '≈ 75 % hielo y 25 % roca', atm: 'Exosfera tenue de oxígeno y dióxido de carbono',
feats: ['Superficie antigua cubierta de cráteres.', 'Grietas brillantes en su hemisferio posterior.'],
facts: ['Cassini detectó su exosfera en 2010.'],
}, src: SRC_MOON,
},
{
id: 'titan', name: 'Titán', aka: ['titan'], type: 'moon', parent: 'saturno',
R: 2574.7, mass: 1.3452e23, locked: true, color: '#e0a65a',
orbit: { t: 'moon', a: 1221870, P: 15.945421, i: 0.348, node: 0, L0: 160 },
vis: { style: 'rock', feat: 'titan', c1: '#d6943f', c2: '#a8692a', c3: '#e8b56a', crater: 0.0, seed: 331.5, atm: { color: [1.0, 0.62, 0.28], h: 0.09, k: 1.6 } },
info: {
desc: 'La única luna con una atmósfera densa, y el único lugar fuera de la Tierra con lagos líquidos en superficie.',
age: '≈ 4,500 millones de años', temp: '≈ −179 °C', grav: 1.352, moons: 0,
comp: 'Roca y hielo, con un océano interno de agua', atm: '≈ 95 % N₂ y ≈ 5 % CH₄; presión 1.5 veces la terrestre',
feats: ['Neblina anaranjada de compuestos orgánicos.', 'Lagos y mares de metano y etano en las regiones polares.'],
facts: ['La sonda Huygens aterrizó en su superficie en 2005.', 'La misión Dragonfly de la NASA enviará un dron a explorarla.'],
}, src: SRC_MOON,
},
{
id: 'japeto', name: 'Jápeto', aka: ['iapetus', 'japeto'], type: 'moon', parent: 'saturno',
R: 734.5, mass: 1.8056e21, locked: true, color: '#a49c90',
orbit: { t: 'moon', a: 3560820, P: 79.3215, i: 15.47, node: 0, L0: 40 },
vis: { style: 'rock', feat: 'japeto', c1: '#d8d3c8', c2: '#3a2a20', c3: '#efece6', crater: 0.8, seed: 341.2 },
info: {
desc: 'La luna de dos caras: un hemisferio oscuro como el carbón y otro brillante como la nieve.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Principalmente hielo, con poca roca', atm: 'Ninguna',
feats: ['El hemisferio delantero está cubierto de material oscuro.', 'Una cresta ecuatorial de hasta unos 20 km de altura.'],
facts: ['Su órbita está inclinada unos 15° respecto al plano de los anillos.', 'Giovanni Cassini notó en 1671 que sólo era visible en un lado de su órbita.'],
}, src: SRC_MOON,
},
{
id: 'miranda', name: 'Miranda', aka: ['miranda'], type: 'moon', parent: 'urano',
R: 235.8, mass: 6.4e19, locked: true, color: '#b9b6b1',
orbit: { t: 'moon', a: 129390, P: 1.413479, i: 4.232, node: 0, L0: 70 },
vis: { style: 'rock', feat: 'miranda', c1: '#a9a6a1', c2: '#6f6c68', c3: '#d2cfca', crater: 0.6, seed: 351.7 },
info: {
desc: 'Una pequeña luna con un mosaico de terrenos de aspecto caótico.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Hielo de agua y roca', atm: 'Ninguna',
feats: ['Verona Rupes, un acantilado de unos 20 km de altura.', 'Grandes regiones con surcos llamadas coronas.'],
facts: ['Gerard Kuiper la descubrió en 1948.'],
}, src: SRC_MOON,
},
{
id: 'ariel', name: 'Ariel', aka: ['ariel'], type: 'moon', parent: 'urano',
R: 578.9, mass: 1.251e21, locked: true, color: '#c8c6c2',
orbit: { t: 'moon', a: 191020, P: 2.520379, i: 0.26, node: 0, L0: 150 },
vis: { style: 'rock', feat: 'icy', c1: '#bfbcb7', c2: '#85827d', c3: '#e0dedb', crater: 0.6, seed: 361.4 },
info: {
desc: 'La más brillante de las grandes lunas de Urano, con cañones y llanuras jóvenes.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Hielo y roca a partes similares', atm: 'Ninguna',
feats: ['La superficie más joven entre las grandes lunas de Urano.', 'Grandes cañones y llanuras lisas.'],
facts: ['El telescopio James Webb detectó hielo de CO₂ en su superficie.'],
}, src: SRC_MOON,
},
{
id: 'umbriel', name: 'Umbriel', aka: ['umbriel'], type: 'moon', parent: 'urano',
R: 584.7, mass: 1.275e21, locked: true, color: '#7d7a76',
orbit: { t: 'moon', a: 266000, P: 4.144177, i: 0.128, node: 0, L0: 290 },
vis: { style: 'rock', feat: 'generic', c1: '#6a6764', c2: '#474543', c3: '#a7a39e', crater: 0.9, seed: 371.1 },
info: {
desc: 'La más oscura de las grandes lunas de Urano.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Hielo y roca', atm: 'Ninguna',
feats: ['Superficie muy antigua y craterizada.', 'Un anillo brillante en el cráter Wunda.'],
facts: ['William Lassell la descubrió en 1851.'],
}, src: SRC_MOON,
},
{
id: 'titania', name: 'Titania', aka: ['titania'], type: 'moon', parent: 'urano',
R: 788.4, mass: 3.4e21, locked: true, color: '#b2aea8',
orbit: { t: 'moon', a: 435910, P: 8.705872, i: 0.34, node: 0, L0: 20 },
vis: { style: 'rock', feat: 'icy', c1: '#a8a49e', c2: '#706c67', c3: '#cdc9c3', crater: 0.7, seed: 381.8 },
info: {
desc: 'La luna más grande de Urano, marcada por grandes fallas.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Hielo y roca', atm: 'Ninguna',
feats: ['Grandes fallas y cañones como Messina Chasmata.', 'Mezcla de terrenos antiguos y renovados.'],
facts: ['William Herschel la descubrió en 1787.'],
}, src: SRC_MOON,
},
{
id: 'oberon', name: 'Oberón', aka: ['oberon'], type: 'moon', parent: 'urano',
R: 761.4, mass: 3.08e21, locked: true, color: '#a19a92',
orbit: { t: 'moon', a: 583520, P: 13.463239, i: 0.058, node: 0, L0: 230 },
vis: { style: 'rock', feat: 'generic', c1: '#958d84', c2: '#5e5852', c3: '#b9b2a9', crater: 1.0, seed: 391.4 },
info: {
desc: 'La luna principal más externa de Urano, antigua y muy craterizada.',
age: '≈ 4,500 millones de años', temp: 'No disponible', moons: 0,
comp: 'Hielo y roca', atm: 'Ninguna',
feats: ['Cráteres con fondos oscuros.', 'Una montaña de unos 6 km de altura en el limbo.'],
facts: ['William Herschel la descubrió en 1787, junto con Titania.'],
}, src: SRC_MOON,
},
{
id: 'triton', name: 'Tritón', aka: ['triton'], type: 'moon', parent: 'neptuno',
R: 1353.4, mass: 2.139e22, locked: true, color: '#d8c8c0',
orbit: { t: 'moon', a: 354759, P: 5.876854, i: 156.885, node: 0, L0: 300 },
vis: { style: 'rock', feat: 'triton', c1: '#cdbbb2', c2: '#8f7f78', c3: '#f0e6e0', crater: 0.15, seed: 401.1, atm: { color: [0.6, 0.7, 1.0], h: 0.02, k: 0.2 } },
info: {
desc: 'La gran luna de Neptuno, probablemente un objeto del cinturón de Kuiper capturado.',
age: '≈ 4,500 millones de años', temp: '≈ −235 °C', grav: 0.779, moons: 0,
comp: 'Roca y hielo; superficie de nitrógeno congelado', atm: 'Tenue, de nitrógeno',
feats: ['Única luna grande con órbita retrógrada.', 'Géiseres de nitrógeno observados por Voyager 2.'],
facts: ['Es una de las superficies más frías medidas en el Sistema Solar.', 'Su órbita decae lentamente; en miles de millones de años podría romperse por las mareas.'],
}, src: SRC_MOON,
},
{
id: 'caronte', name: 'Caronte', aka: ['charon', 'caronte'], type: 'moon', parent: 'pluton',
R: 606, mass: 1.586e21, locked: true, color: '#a8a4a0',
orbit: { t: 'moon', a: 19591, P: 6.387221, i: 0.0, node: 0, L0: 0 },
vis: { style: 'rock', feat: 'caronte', c1: '#9c9893', c2: '#6b4a3a', c3: '#bdb9b4', crater: 0.6, seed: 411.6 },
info: {
desc: 'La mayor luna de Plutón, tan grande que ambos forman casi un sistema doble.',
age: '≈ 4,500 millones de años', temp: '≈ −220 °C', grav: 0.288, moons: 0,
comp: 'Hielo de agua y roca', atm: 'Ninguna',
feats: ['Región polar rojiza (Mordor Macula) teñida por material escapado de Plutón.', 'Cañones de hasta unos 9 km de profundidad.'],
facts: ['Plutón y Caronte giran alrededor de un punto situado fuera de Plutón.', 'Ambos se muestran siempre la misma cara.'],
}, src: SRC_MOON + ' · New Horizons (NASA)',
},
];
const BODY = {};
BODIES.forEach(b => { BODY[b.id] = b; b.children = []; });
BODIES.forEach(b => { if (b.parent && BODY[b.parent]) BODY[b.parent].children.push(b); });
const TOUR = [
{ id: 'sol', txt: 'El Sol contiene el 99.86 % de la masa del Sistema Solar.' },
{ id: 'mercurio', txt: 'Mercurio gira tres veces por cada dos vueltas al Sol.' },
{ id: 'venus', txt: 'En Venus la presión en superficie es unas 92 veces la terrestre.' },
{ id: 'tierra', txt: 'La Tierra: océanos, nubes y el único lugar conocido con vida.' },
{ id: 'luna', txt: 'La Luna siempre muestra la misma cara a la Tierra.' },
{ id: 'marte', txt: 'Marte alberga Olympus Mons, el mayor volcán conocido.' },
{ id: 'ceres', txt: 'Ceres, el planeta enano del cinturón de asteroides.' },
{ id: 'jupiter', txt: 'Júpiter: más del doble de masivo que el resto de planetas juntos.' },
{ id: 'io', txt: 'Ío es el cuerpo más volcánico del Sistema Solar.' },
{ id: 'saturno', txt: 'Los anillos de Saturno son de hielo y apenas tienen decenas de metros de espesor.' },
{ id: 'titan', txt: 'Titán tiene lagos de metano bajo una atmósfera densa.' },
{ id: 'urano', txt: 'Urano gira con el eje inclinado casi 98°.' },
{ id: 'neptuno', txt: 'En Neptuno soplan vientos de más de 2,000 km/h.' },
{ id: 'triton', txt: 'Tritón orbita a Neptuno en sentido contrario a su rotación.' },
{ id: 'pluton', txt: 'Plutón y su corazón de hielo de nitrógeno.' },
{ id: 'halley', txt: 'El cometa Halley regresará al Sistema Solar interior hacia 2061.' },
];
const SRC_CRAFT = 'NASA, ESA, JAXA, CSA, CMSA (páginas oficiales de cada misión)';
const CRAFTS = [
{
id: 'jwst', name: 'Telescopio espacial James Webb', short: 'James Webb', aka: ['jwst', 'webb', 'james webb', 'telescopio webb'], type: 'craft', parent: 'tierra',
sub: 'Telescopio espacial infrarrojo', R: 0.0106, color: '#f2c879', model: 'jwst', validFrom: '2022-01-24',
orbit: { t: 'lpoint', L: 2, dist: 1.5e6, halo: [520000, 800000], P: 180, phase: 40 },
craft: { agency: 'NASA, ESA y CSA', launch: '2021-12-25', launchTxt: '25 de diciembre de 2021 (Ariane 5, Kurú)', orbitTxt: 'Órbita halo alrededor del punto L2 Sol-Tierra, a unos 1.5 millones de km de la Tierra', periodTxt: '≈ 6 meses por vuelta al halo', dims: 'Parasol de 21.2 × 14.2 m; espejo primario de 6.5 m', mass: '≈ 6,200 kg', goal: 'Observar en infrarrojo las primeras galaxias, la formación de estrellas y planetas, y las atmósferas de exoplanetas', status: 'En operación' },
info: {
desc: 'El mayor telescopio espacial jamás lanzado. Observa en infrarrojo desde el punto L2, protegido del calor del Sol por un parasol del tamaño de una cancha de tenis.',
feats: ['Espejo primario de 18 segmentos hexagonales de berilio recubiertos de oro.', 'Parasol de 5 capas que mantiene los instrumentos por debajo de unos −223 °C.', 'Desde L2, el Sol, la Tierra y la Luna quedan siempre del mismo lado, tras el parasol.'],
facts: ['Sus primeras imágenes científicas se publicaron el 12 de julio de 2022.', 'El lanzamiento fue tan preciso que ahorró combustible para más de 20 años de misión.'],
}, src: 'NASA/ESA/CSA (Webb)',
},
{
id: 'euclid', name: 'Euclid', short: 'Euclid', aka: ['euclid'], type: 'craft', parent: 'tierra',
sub: 'Telescopio espacial de cartografía cósmica', R: 0.0024, color: '#cfd8e6', model: 'euclid', validFrom: '2023-08-01',
orbit: { t: 'lpoint', L: 2, dist: 1.5e6, halo: [380000, 640000], P: 180, phase: 200 },
craft: { agency: 'ESA (con participación de la NASA)', launch: '2023-07-01', launchTxt: '1 de julio de 2023 (Falcon 9)', orbitTxt: 'Órbita halo alrededor del punto L2 Sol-Tierra', periodTxt: '≈ 6 meses por vuelta al halo', dims: '≈ 4.7 × 3.7 m; telescopio de 1.2 m', mass: '≈ 2,000 kg', goal: 'Cartografiar miles de millones de galaxias para estudiar la materia oscura y la energía oscura', status: 'En operación' },
info: {
desc: 'Un telescopio que está trazando el mapa 3D más grande del universo oscuro, observando más de un tercio del cielo.',
feats: ['Combina una cámara de luz visible y un espectrofotómetro infrarrojo.', 'Mide la forma de miles de millones de galaxias para detectar la materia oscura por su efecto gravitatorio.'],
facts: ['Comparte la región de L2 con el telescopio James Webb.', 'La ESA publicó su primer gran conjunto de datos en marzo de 2025.'],
}, src: 'ESA (Euclid)',
},
{
id: 'soho', name: 'SOHO', short: 'SOHO', aka: ['soho', 'observatorio solar'], type: 'craft', parent: 'tierra',
sub: 'Observatorio solar', R: 0.00475, color: '#e3b65a', model: 'soho', validFrom: '1996-02-14',
orbit: { t: 'lpoint', L: 1, dist: 1.5e6, halo: [200000, 650000], P: 178, phase: 110 },
craft: { agency: 'ESA y NASA', launch: '1995-12-02', launchTxt: '2 de diciembre de 1995 (Atlas IIAS)', orbitTxt: 'Órbita halo alrededor del punto L1 Sol-Tierra, a unos 1.5 millones de km hacia el Sol', periodTxt: '≈ 6 meses por vuelta al halo', dims: 'Cuerpo de 4.3 × 2.7 × 3.7 m; 9.5 m con paneles solares', mass: '≈ 1,850 kg', goal: 'Estudiar el interior, la atmósfera y el viento del Sol', status: 'En operación' },
info: {
desc: 'Un observatorio que vigila el Sol sin interrupción desde el punto L1, entre la Tierra y el Sol.',
feats: ['Su coronógrafo bloquea el disco solar para observar la corona y las eyecciones de masa coronal.', 'Sus datos ayudan a pronosticar el clima espacial.'],
facts: ['Ha descubierto más de 5,000 cometas, la mayoría rozadores del Sol.', 'Se perdió el contacto en 1998 y fue recuperado meses después.'],
}, src: 'ESA/NASA (SOHO)',
},
{
id: 'iss', name: 'Estación Espacial Internacional', short: 'EEI', aka: ['iss', 'eei', 'estacion espacial internacional', 'estacion espacial'], type: 'craft', parent: 'tierra',
sub: 'Estación espacial tripulada', R: 0.0545, color: '#e9edf3', model: 'iss', validFrom: '1998-11-20', validTo: '2031-01-01',
orbit: { t: 'leo', a: 6371 + 415, P: 0.0644, i: 51.64, node: 120, L0: 30 },
craft: { agency: 'NASA, Roscosmos, ESA, JAXA y CSA', launch: '1998-11-20', launchTxt: '20 de noviembre de 1998 (primer módulo, Zarya)', orbitTxt: 'Órbita terrestre baja a unos 400–420 km de altitud, inclinada 51.6°', periodTxt: '≈ 93 minutos', dims: '≈ 109 × 73 m', mass: '≈ 420,000 kg', goal: 'Laboratorio de investigación en microgravedad', status: 'En operación; su retiro está previsto hacia 2030' },
info: {
desc: 'El objeto más grande construido por la humanidad en el espacio, habitado de forma continua desde noviembre de 2000.',
feats: ['Viaja a unos 28,000 km/h y da unas 16 vueltas a la Tierra cada día.', 'Sus paneles solares cubren una superficie mayor que media cancha de fútbol americano.'],
facts: ['Se ensambló en órbita a lo largo de más de una década, con decenas de lanzamientos.', 'Se ve a simple vista como un punto brillante que cruza el cielo en pocos minutos.'],
}, src: 'NASA (ISS)',
},
{
id: 'tiangong', name: 'Estación espacial Tiangong', short: 'Tiangong', aka: ['tiangong', 'estacion china'], type: 'craft', parent: 'tierra',
sub: 'Estación espacial tripulada', R: 0.02, color: '#e9edf3', model: 'tiangong', validFrom: '2021-04-29',
orbit: { t: 'leo', a: 6371 + 390, P: 0.0640, i: 41.5, node: 300, L0: 200 },
craft: { agency: 'Agencia Espacial Tripulada de China (CMSA)', launch: '2021-04-29', launchTxt: '29 de abril de 2021 (módulo central Tianhe)', orbitTxt: 'Órbita terrestre baja a unos 380–390 km de altitud, inclinada 41.5°', periodTxt: '≈ 92 minutos', dims: 'Configuración en T de tres módulos', mass: '≈ 66,000 kg (sin naves acopladas)', goal: 'Laboratorio orbital tripulado', status: 'En operación' },
info: {
desc: 'La estación espacial de China, formada por tres módulos en forma de T.',
feats: ['Módulo central Tianhe y laboratorios Wentian y Mengtian.', 'Completada en noviembre de 2022.'],
facts: ['Su nombre significa "Palacio celestial".'],
}, src: 'CMSA',
},
{
id: 'hubble', name: 'Telescopio espacial Hubble', short: 'Hubble', aka: ['hubble', 'telescopio hubble', 'hst'], type: 'craft', parent: 'tierra',
sub: 'Telescopio espacial', R: 0.0066, color: '#d6dae0', model: 'hubble', validFrom: '1990-04-25',
orbit: { t: 'leo', a: 6371 + 515, P: 0.0660, i: 28.47, node: 30, L0: 110 },
craft: { agency: 'NASA y ESA', launch: '1990-04-24', launchTxt: '24 de abril de 1990 (transbordador Discovery)', orbitTxt: 'Órbita terrestre baja a unos 500 km de altitud, inclinada 28.5°', periodTxt: '≈ 95 minutos', dims: '13.2 m de largo; espejo de 2.4 m', mass: '≈ 12,000 kg', goal: 'Observar en luz ultravioleta, visible e infrarroja cercana', status: 'En operación' },
info: {
desc: 'El telescopio espacial que transformó la astronomía moderna, en órbita desde 1990.',
feats: ['Fue diseñado para recibir mantenimiento de astronautas en órbita.', 'Ayudó a medir la expansión acelerada del universo.'],
facts: ['Cinco misiones de servicio, entre 1993 y 2009, lo repararon y mejoraron.', 'La primera de ellas corrigió el defecto de su espejo principal.'],
}, src: 'NASA/ESA (Hubble)',
},
{
id: 'parker', name: 'Sonda solar Parker', short: 'Parker', aka: ['parker', 'parker solar probe', 'sonda parker'], type: 'craft', parent: 'sol',
sub: 'Sonda solar', R: 0.0015, color: '#f0f0ea', model: 'parker', validFrom: '2024-11-07',
orbit: { t: 'kep', a: 0.38725, e: 0.8815, i: 3.4, node: 76.7, peri: 30, tp: '2024-12-24' },
craft: { agency: 'NASA', launch: '2018-08-12', launchTxt: '12 de agosto de 2018 (Delta IV Heavy)', orbitTxt: 'Órbita alrededor del Sol entre unos 6.9 millones de km y la órbita de Venus', periodTxt: '≈ 88 días', dims: 'Escudo térmico de 2.4 m de diámetro', mass: '≈ 685 kg', goal: 'Atravesar la corona solar para estudiar el origen del viento solar', status: 'En órbita solar' },
info: {
desc: 'La nave que más se ha acercado al Sol y el objeto más rápido construido por el ser humano.',
feats: ['Su escudo de compuesto de carbono resiste unos 1,400 °C mientras el interior se mantiene a temperatura ambiente.', 'En cada perihelio atraviesa la corona del Sol.'],
facts: ['El 24 de diciembre de 2024 pasó a unos 6.1 millones de km de la superficie solar.', 'En ese paso alcanzó unos 692,000 km/h.'],
}, src: 'NASA (Parker Solar Probe)',
},
{
id: 'voyager1', name: 'Voyager 1', short: 'Voyager 1', aka: ['voyager 1', 'voyager'], type: 'craft', parent: 'sol',
sub: 'Sonda interestelar', R: 0.0065, color: '#dfe6f0', model: 'voyager', validFrom: '1981-01-01',
orbit: { t: 'drift', ra: 258.3, dec: 12.0, r0: 163.82, t0: '2024-06-15', rate: 3.58 },
craft: { agency: 'NASA/JPL', launch: '1977-09-05', launchTxt: '5 de septiembre de 1977 (Titan IIIE)', orbitTxt: 'Trayectoria de escape del Sistema Solar, en el espacio interestelar desde 2012', periodTxt: 'No aplica: no regresará', dims: 'Antena de 3.7 m; brazo magnetométrico de 13 m', mass: '≈ 825 kg al lanzamiento', goal: 'Explorar Júpiter y Saturno; hoy estudia el medio interestelar', status: 'En operación con instrumentos reducidos para ahorrar energía' },
info: {
desc: 'El objeto más lejano construido por la humanidad, viajando por el espacio interestelar.',
feats: ['Cruzó la heliopausa el 25 de agosto de 2012.', 'Lleva el Disco de Oro con sonidos e imágenes de la Tierra.'],
facts: ['La NASA estima que alcanzará un día luz de distancia de la Tierra el 18 de noviembre de 2026.', 'En 1990 tomó la fotografía de la Tierra conocida como "Un punto azul pálido".'],
}, src: 'NASA/JPL (Voyager)',
},
{
id: 'voyager2', name: 'Voyager 2', short: 'Voyager 2', aka: ['voyager 2'], type: 'craft', parent: 'sol',
sub: 'Sonda interestelar', R: 0.0065, color: '#dfe6f0', model: 'voyager', validFrom: '1990-01-01',
orbit: { t: 'drift', ra: 300.0, dec: -59.0, r0: 136.81, t0: '2024-06-15', rate: 3.17 },
craft: { agency: 'NASA/JPL', launch: '1977-08-20', launchTxt: '20 de agosto de 1977 (Titan IIIE)', orbitTxt: 'Trayectoria de escape del Sistema Solar, en el espacio interestelar desde 2018', periodTxt: 'No aplica: no regresará', dims: 'Antena de 3.7 m; brazo magnetométrico de 13 m', mass: '≈ 825 kg al lanzamiento', goal: 'Explorar los cuatro planetas gigantes; hoy estudia el medio interestelar', status: 'En operación con instrumentos reducidos para ahorrar energía' },
info: {
desc: 'La única nave que ha visitado Urano y Neptuno, hoy en el espacio interestelar.',
feats: ['Sobrevoló Júpiter (1979), Saturno (1981), Urano (1986) y Neptuno (1989).', 'Cruzó la heliopausa el 5 de noviembre de 2018.'],
facts: ['Se lanzó 16 días antes que su gemela, la Voyager 1.'],
}, src: 'NASA/JPL (Voyager)',
},
{
id: 'newhorizons', name: 'New Horizons', short: 'New Horizons', aka: ['new horizons', 'nuevos horizontes'], type: 'craft', parent: 'sol',
sub: 'Sonda del cinturón de Kuiper', R: 0.0014, color: '#e8c87a', model: 'newhorizons', validFrom: '2016-01-01',
orbit: { t: 'drift', ra: 290.5, dec: -20.5, r0: 58.6, t0: '2024-01-01', rate: 2.87 },
craft: { agency: 'NASA', launch: '2006-01-19', launchTxt: '19 de enero de 2006 (Atlas V)', orbitTxt: 'Trayectoria de escape a través del cinturón de Kuiper', periodTxt: 'No aplica: no regresará', dims: '≈ 2.7 × 2.1 m; antena de 2.1 m', mass: '≈ 478 kg al lanzamiento', goal: 'Explorar Plutón y objetos del cinturón de Kuiper', status: 'En misión extendida en el cinturón de Kuiper' },
info: {
desc: 'La primera nave en explorar Plutón y el objeto más lejano jamás visitado, Arrokoth.',
feats: ['Sobrevoló Plutón el 14 de julio de 2015.', 'Sobrevoló Arrokoth el 1 de enero de 2019.'],
facts: ['Lleva parte de las cenizas de Clyde Tombaugh, el descubridor de Plutón.'],
}, src: 'NASA/JHUAPL (New Horizons)',
},
];
CRAFTS.forEach(c => {
c.vis = { style: 'craft' };
BODIES.push(c); BODY[c.id] = c; c.children = [];
BODY[c.parent].children.push(c);
});
TYPE_LABEL.craft = 'Nave espacial';
TOUR.splice(5, 0, { id: 'jwst', txt: 'El James Webb observa en infrarrojo desde L2, a 1.5 millones de km de la Tierra.' });
TOUR.splice(TOUR.length - 1, 0, { id: 'voyager1', txt: 'Voyager 1, lanzada en 1977, es el objeto humano más lejano: viaja por el espacio interestelar.' });
const JWST_PACK = 'assets/models/satellites/jwst.bin';
const JWST_PALETTE = [[0.554, 0.554, 0.554, 0.745], [0.623, 0.427, 0.574, 0.75], [0.827, 0.828, 0.805, 0.262], [0.462, 0.462, 0.462, 0.35], [0.367, 0.367, 0.367, 0.63], [0.122, 0.122, 0.122, 0.262], [0.392, 0.392, 0.392, 0.63], [0.778, 0.74, 0.574, 0.262], [0.696, 0.697, 0.678, 0.262], [0.827, 0.828, 0.805, 0.262], [0.656, 0.45, 0.605, 0.262], [0.623, 0.427, 0.574, 0.75], [0.087, 0.08, 0.1, 0.262], [0.2, 0.2, 0.2, 0.262], [0.15, 0.105, 0.12, 0.262], [0.765, 0.659, 0.264, 0.262], [0.85, 0.85, 0.85, 0.3], [0.15, 0.105, 0.12, 0.262], [0.92, 0.66, 0.2, 1.5], [0.401, 0.403, 0.406, 0.262], [0.652, 0.652, 0.652, 0.262], [0.92, 0.66, 0.2, 1.5], [0.42, 0.36, 0.6, 0.8], [0.746, 0.746, 0.746, 0.262], [0.479, 0.479, 0.479, 0.262], [0.246, 0.246, 0.246, 0.262], [0.78, 0.55, 0.17, 0.6], [0.46, 0.36, 0.5, 0.3], [0.5, 0.47, 0.55, 0.45], [0.48, 0.45, 0.53, 0.45], [0.47, 0.44, 0.52, 0.45], [0.5, 0.47, 0.56, 0.45], [0.15, 0.138, 0.075, 0.262], [0.09, 0.11, 0.3, 0.7], [0.623, 0.427, 0.574, 0.75], [0.784, 0.784, 0.784, 0.262]];
const MODEL_CREDITS = {
jwst: { title: 'JWST (james webb space telescope)', author: 'Paul', url: 'https://sketchfab.com/3d-models/jwst-james-webb-space-telescope-6c92c08a672640afb58ee44d248fd0fe', license: 'CC BY 4.0', licenseUrl: 'https://creativecommons.org/licenses/by/4.0/', note: 'adaptado: escala, materiales y formato' },
};
MODEL_CREDITS.hubble = { title: 'Hubble', author: 'uperesito', url: 'https://sketchfab.com/3d-models/hubble-1b82249a0c3540daa66092a8d73442e9', license: 'CC BY 4.0', licenseUrl: 'https://creativecommons.org/licenses/by/4.0/', note: 'adaptado: escala, simplificación y formato' };
MODEL_CREDITS.voyager = { title: 'Voyager 1', author: 'Sidharth1', url: 'https://sketchfab.com/3d-models/voyager-1-f626e0b2d37b46dab4c9433306720c4b', license: 'CC BY 4.0', licenseUrl: 'https://creativecommons.org/licenses/by/4.0/', note: 'adaptado: escala, orientación y formato' };
MODEL_CREDITS.iss = { title: 'Estación Espacial Internacional', author: null, url: null, license: null, licenseUrl: null, note: 'modelo aportado por el usuario; el archivo no indica autor ni licencia' };
BODY.jwst.credit = MODEL_CREDITS.jwst;
MODEL_CREDITS.tiangong = { title: 'Estación espacial Tiangong', author: null, url: null, license: null, licenseUrl: null, note: 'modelo aportado por el creador del proyecto y completado para SOLARIS: paneles solares del laboratorio inferior, normales de superficie y escala real' };
BODY.hubble.credit = MODEL_CREDITS.hubble; BODY.iss.credit = MODEL_CREDITS.iss; BODY.tiangong.credit = MODEL_CREDITS.tiangong;
BODY.voyager1.credit = MODEL_CREDITS.voyager; BODY.voyager2.credit = MODEL_CREDITS.voyager;
const PLANET_TEX = {"sol": {"credit": {"title": "Sun", "author": "SebastianSosnowski", "url": "https://sketchfab.com/3d-models/sun-9ef1c68fbb944147bcfcc891d3912645"}, "map": "assets/textures/planets/sol.webp"}, "mercurio": {"credit": {"title": "Mercurio v1.1", "author": "uperesito", "url": "https://sketchfab.com/3d-models/mercurio-v11-06baa7da7c9743dc9ef0111e1db0e8eb", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}, "prevCredit": {"title": "Mercury", "author": "Akshat", "url": "https://sketchfab.com/3d-models/mercury-32347fa4ec1a4987b71f461a401d91c4"}, "map": "assets/textures/planets/mercurio.webp"}, "venus": {"credit": {"title": "Venus", "author": "Akshat", "url": "https://sketchfab.com/3d-models/venus-d497ce25553447f3b7b679110c85cfa1"}, "map": "assets/textures/planets/venus.webp"}, "tierra": {"credit": {"title": "Earth", "author": "Akshat", "url": "https://sketchfab.com/3d-models/earth-41fc80d85dfd480281f21b74b2de2faa"}, "map": "assets/textures/planets/tierra.webp"}, "luna": {"nScale": 0.43, "credit": {"title": "The Moon", "author": "SebastianSosnowski", "url": "https://sketchfab.com/3d-models/the-moon-9916fcec59f04b07b3e8d7f077dc3ded"}, "map": "assets/textures/planets/luna.webp", "nmap": "assets/textures/planets/luna-relieve.jpg"}, "marte": {"credit": {"title": "Mars", "author": "Akshat", "url": "https://sketchfab.com/3d-models/mars-9c7bbc64d8c74acfa9ec344c0fc10e1a"}, "map": "assets/textures/planets/marte.webp"}, "jupiter": {"credit": {"title": "Realistic Jupiter", "author": "Shady Tex", "url": "https://sketchfab.com/3d-models/realistic-jupiter-993ba62a539e4c308e9e3137df454ed6", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}, "prevCredit": {"title": "Jupiter", "author": "Akshat", "url": "https://sketchfab.com/3d-models/jupiter-d252c96ae3de48d7968b1206522ba9f5"}, "map": "assets/textures/planets/jupiter.webp"}, "saturno": {"credit": {"title": "Saturn", "author": "Akshat", "url": "https://sketchfab.com/3d-models/saturn-8fb67d3defd74aaa880df3a08317e641"}, "map": "assets/textures/planets/saturno.webp", "ring": {"uri": "assets/textures/planets/saturno-anillos.png", "rin": 1.239, "rout": 2.352}}, "urano": {"credit": {"title": "Uranus", "author": "Akshat", "url": "https://sketchfab.com/3d-models/uranus-4d2f0c3674904472ac413fdabbf491d7"}, "map": "assets/textures/planets/urano.webp"}, "neptuno": {"credit": {"title": "Neptune", "author": "Akshat", "url": "https://sketchfab.com/3d-models/neptune-947a405a0a4348f9a49ff4bd3ed3cc4b"}, "map": "assets/textures/planets/neptuno.webp"}, "pluton": {"credit": {"title": "Pluto", "author": "Akshat", "url": "https://sketchfab.com/3d-models/pluto-a30ff3a5f4f3477d87fc534e0d1df7e2"}, "map": "assets/textures/planets/pluton.webp"}, "ganimedes": {"nScale": 0.3, "credit": {"title": "Ganímedes v1.1", "author": "uperesito", "url": "https://sketchfab.com/3d-models/ganimedes-v11-e991ac7326984467b172fa1723edff93"}, "map": "assets/textures/planets/ganimedes.webp", "nmap": "assets/textures/planets/ganimedes-relieve.jpg"}, "titan": {"credit": {"title": "Titan", "author": "uperesito", "url": "https://sketchfab.com/3d-models/titan-f5f5989b3ebf4d7a96e52a49ab6966ef"}, "map": "assets/textures/planets/titan.webp"}, "io": {"map": "assets/textures/moons/io.webp", "credit": {"title": "Io v1.1", "author": "uperesito", "url": "https://sketchfab.com/3d-models/io-v11-24eff42e8ac2429a945fddd722f4ef5c", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}, "nmap": "assets/textures/moons/io-relieve.jpg", "nScale": 0.3}, "calisto": {"map": "assets/textures/moons/calisto.webp", "credit": {"title": "Calisto v1.1", "author": "uperesito", "url": "https://sketchfab.com/3d-models/calisto-v11-ca28b6b42fcb4e589ca5c830bc4e6427", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}, "nmap": "assets/textures/moons/calisto-relieve.jpg", "nScale": 0.3}, "rea": {"map": "assets/textures/moons/rea.webp", "credit": {"title": "Rea", "author": "uperesito", "url": "https://sketchfab.com/3d-models/rea-ea0ead107be644c999e6f4f21d9e4970", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}}, "oberon": {"map": "assets/textures/moons/oberon.webp", "credit": {"title": "Oberon", "author": "uperesito", "url": "https://sketchfab.com/3d-models/oberon-97a4c117709c48e68a044254c9820159", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}}, "ariel": {"map": "assets/textures/moons/ariel.webp", "credit": {"title": "Ariel", "author": "uperesito", "url": "https://sketchfab.com/3d-models/ariel-26b522635a764e68b4971a1879998ede", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}}, "miranda": {"map": "assets/textures/moons/miranda.webp", "credit": {"title": "Miranda", "author": "uperesito", "url": "https://sketchfab.com/3d-models/miranda-40106d3936544c6b89209ff31a82869a", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: orientación normalizada, tamaño y formato"}}, "fobos": {"map": "assets/textures/moons/fobos.webp", "credit": {"title": "Fobos", "author": "uperesito", "url": "https://sketchfab.com/3d-models/fobos-5226064b82274ce7aae865dc1cbae9ac", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "forma y textura; adaptado: escala, simplificación de la geometría y formato"}}, "deimos": {"map": "assets/textures/moons/deimos.webp", "credit": {"title": "Deimos", "author": "uperesito", "url": "https://sketchfab.com/3d-models/deimos-db9c281d0a7b4452aaac26eb1d3738ad", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "forma y textura; adaptado: escala, simplificación de la geometría y formato"}}, "ceres": {"map": "assets/textures/dwarfs/ceres.webp", "credit": {"title": "Ceres", "author": "Mieke Roth", "url": "https://sketchfab.com/3d-models/ceres-c3ee7020e07e47e9b3d74585ad612bb9", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: proyección y orientación normalizadas, tamaño y formato"}}, "haumea": {"map": "assets/textures/dwarfs/haumea.webp", "credit": {"title": "Haumea", "author": "ПГГПУ", "url": "https://sketchfab.com/3d-models/haumea-3f93a3b07cce45cd809f0521d8ef75cd", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "forma (elipsoide de 2322 × 1704 × 1138 km) y textura; adaptado: remallado a mayor resolución, orientación normalizada y formato"}, "nmap": "assets/textures/dwarfs/haumea-relieve.jpg", "nScale": 0.3}, "makemake": {"map": "assets/textures/dwarfs/makemake.webp", "credit": {"title": "Makemake (Got Balls Remastered)", "author": "stanimation3d", "url": "https://sketchfab.com/3d-models/makemake-got-balls-remastered-f9e5ceacb4bd4aeca132193fd7e134f0", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: proyección y orientación normalizadas, tamaño y formato"}}, "sedna": {"map": "assets/textures/dwarfs/sedna.webp", "credit": {"title": "Sedna (Got Balls Remastered)", "author": "stanimation3d", "url": "https://sketchfab.com/3d-models/sedna-got-balls-remastered-049c9b0628b44735b8a1833dd58af8dc", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura adaptada: proyección y orientación normalizadas, tamaño y formato"}}, "quaoar": {"map": "assets/textures/dwarfs/quaoar.webp", "credit": {"title": "Quaoar (Got Balls Remastered)", "author": "stanimation3d", "url": "https://sketchfab.com/3d-models/quaoar-got-balls-remastered-fa9a458ef8c2434794ee54066db01be9", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "textura y anillo; adaptado: proyección y orientación normalizadas, anillo a partir del perfil del disco del modelo, y formato"}, "rings": {"in": 6.01, "out": 7.19, "type": "faint", "alpha": 0.486}}};
Object.keys(PLANET_TEX).forEach(id => { const c = PLANET_TEX[id].credit; BODY[id].credit = { title: c.title, author: c.author, url: c.url, license: 'CC BY 4.0', licenseUrl: 'https://creativecommons.org/licenses/by/4.0/', note: 'textura adaptada: orientación normalizada y tamaño' }; BODY[id].texModel = true; });
BODY.venus.info.feats.push('La textura muestra la superficie según los mapas de radar: a simple vista, las nubes la ocultan por completo.');
BODY.titan.info.feats.push('Su superficie se representa a partir de observaciones en infrarrojo y radar; en luz visible la neblina anaranjada la oculta casi por completo.');
const MODEL_PACKS = {"voyager": {"mats": [{"name": "material_1", "col": [1, 1, 1], "spec": 0.251, "tex": 0}, {"name": "material_2", "col": [1.0, 0.827, 0.318], "spec": 0.608, "tex": -1}, {"name": "material_3", "col": [1, 1, 1], "spec": 0.251, "tex": 1}, {"name": "material_4", "col": [1, 1, 1], "spec": 0.251, "tex": 2}], "credit": {"title": "Voyager 1", "author": "Sidharth1", "url": "https://sketchfab.com/3d-models/voyager-1-f626e0b2d37b46dab4c9433306720c4b"}, "pack": "assets/models/satellites/voyager.bin", "tex": ["assets/textures/spacecraft/voyager-0.avif", "assets/textures/spacecraft/voyager-1.avif", "assets/textures/spacecraft/voyager-2.avif"], "texAlt": [null, null, null]}, "hubble": {"mats": [{"name": "Material9", "col": [1, 1, 1], "spec": 0.42, "tex": 0}, {"name": "Material", "col": [1, 1, 1], "spec": 0.42, "tex": 1}, {"name": "Material10", "col": [1, 1, 1], "spec": 0.42, "tex": 2}, {"name": "Material3", "col": [1, 1, 1], "spec": 0.42, "tex": 3}, {"name": "Material11", "col": [1, 1, 1], "spec": 0.42, "tex": 4}, {"name": "Material12", "col": [1, 1, 1], "spec": 0.42, "tex": 5}, {"name": "Material4", "col": [1, 1, 1], "spec": 0.42, "tex": 6}, {"name": "Material6", "col": [1, 1, 1], "spec": 0.42, "tex": 6}, {"name": "Material8", "col": [0.542, 0.542, 0.542], "spec": 0.42, "tex": -1}, {"name": "Material7", "col": [1, 1, 1], "spec": 0.225, "tex": 7}, {"name": "scblack", "col": [1, 1, 1], "spec": 0.42, "tex": 8}, {"name": "Material2", "col": [1, 1, 1], "spec": 0.42, "tex": 9}, {"name": "al_irradite", "col": [1, 1, 1], "spec": 0.42, "tex": 10}, {"name": "logos", "col": [1, 1, 1], "spec": 0.42, "tex": 11}, {"name": "GLOBAL03", "col": [0.3, 0.3, 0.3], "spec": 0.42, "tex": -1}, {"name": "scblack1", "col": [1, 1, 1], "spec": 0.42, "tex": 8}, {"name": "GLOBAL04", "col": [1, 1, 1], "spec": 0.42, "tex": 9}], "credit": {"title": "Hubble", "author": "uperesito", "url": "https://sketchfab.com/3d-models/hubble-1b82249a0c3540daa66092a8d73442e9"}, "pack": "assets/models/satellites/hubble.bin", "tex": ["assets/textures/spacecraft/hubble-0.jpg", "assets/textures/spacecraft/hubble-1.avif", "assets/textures/spacecraft/hubble-2.jpg", "assets/textures/spacecraft/hubble-3.jpg", "assets/textures/spacecraft/hubble-4.avif", "assets/textures/spacecraft/hubble-5.avif", "assets/textures/spacecraft/hubble-6.avif", "assets/textures/spacecraft/hubble-7.jpg", "assets/textures/spacecraft/hubble-8.jpg", "assets/textures/spacecraft/hubble-9.jpg", "assets/textures/spacecraft/hubble-10.jpg", "assets/textures/spacecraft/hubble-11.jpg"], "texAlt": [null, null, null, null, null, null, null, null, null, null, null, null]}, "iss": {"mats": [{"name": "material_0", "col": [0.85, 0.85, 0.85], "spec": 0.12, "tex": -1}, {"name": "ISS_03_dull", "col": [1, 1, 1], "spec": 0.12, "tex": 0}, {"name": "ISS_AO_06", "col": [1, 1, 1], "spec": 0.398, "tex": 1}, {"name": "ISS_01_dull", "col": [1, 1, 1], "spec": 0.366, "tex": 2}, {"name": "ISS_02_dull", "col": [1, 1, 1], "spec": 0.347, "tex": 3}, {"name": "shiny_panel", "col": [1, 1, 1], "spec": 0.443, "tex": 3}, {"name": "olive *", "col": [0.432, 0.425, 0.389], "spec": 0.645, "tex": -1}, {"name": "ISS_AO_02", "col": [1, 1, 1], "spec": 0.533, "tex": 4}, {"name": "ISS_AO_04", "col": [1, 1, 1], "spec": 0.398, "tex": 5}, {"name": "ISS_AO_01", "col": [1, 1, 1], "spec": 0.398, "tex": 6}, {"name": "ISS_AO_03", "col": [1, 1, 1], "spec": 0.398, "tex": 7}, {"name": "ISS_01_shiny_n", "col": [1, 1, 1], "spec": 0.398, "tex": 2}, {"name": "ISS_01_dark *", "col": [1, 1, 1], "spec": 0.251, "tex": 2}, {"name": "foil_silver", "col": [0.85, 0.85, 0.85], "spec": 0.645, "tex": -1}, {"name": "ISS_AO_05", "col": [1, 1, 1], "spec": 0.398, "tex": 8}, {"name": "ISS_04_dull", "col": [1, 1, 1], "spec": 0.251, "tex": 9}, {"name": "ISS_02_dark *", "col": [1, 1, 1], "spec": 0.251, "tex": 3}, {"name": "ISS_AO_07", "col": [1, 1, 1], "spec": 0.398, "tex": 10}, {"name": "ISS_AO_08", "col": [1, 1, 1], "spec": 0.398, "tex": 11}, {"name": "ISS_03_shiny_n", "col": [1, 1, 1], "spec": 0.398, "tex": 0}, {"name": "ecostress", "col": [1, 1, 1], "spec": 0.186, "tex": 12}, {"name": "plastic black", "col": [0.026, 0.026, 0.026], "spec": 0.251, "tex": -1}, {"name": "ISS_03_dull.002", "col": [1, 1, 1], "spec": 0.186, "tex": 13}, {"name": "ecostress metal", "col": [1, 1, 1], "spec": 0.259, "tex": 12}, {"name": "ecostress_dexter", "col": [1, 1, 1], "spec": 0.251, "tex": 2}, {"name": "white", "col": [0.497, 0.497, 0.497], "spec": 0.159, "tex": -1}, {"name": "base_metal", "col": [0.475, 0.473, 0.439], "spec": 0.12, "tex": -1}, {"name": "ISS_03_dull.001", "col": [1, 1, 1], "spec": 0.365, "tex": 14}], "credit": {"title": "EEI (archivo proporcionado por el usuario)", "author": null, "url": null}, "pack": "assets/models/satellites/iss.bin", "tex": ["assets/textures/spacecraft/iss-0.avif", "assets/textures/spacecraft/iss-1.avif", "assets/textures/spacecraft/iss-2.avif", "assets/textures/spacecraft/iss-3.avif", "assets/textures/spacecraft/iss-4.avif", "assets/textures/spacecraft/iss-5.avif", "assets/textures/spacecraft/iss-6.avif", "assets/textures/spacecraft/iss-7.avif", "assets/textures/spacecraft/iss-8.avif", "assets/textures/spacecraft/iss-9.avif", "assets/textures/spacecraft/iss-10.avif", "assets/textures/spacecraft/iss-11.avif", "assets/textures/spacecraft/iss-12.jpg", "assets/textures/spacecraft/iss-13.avif", "assets/textures/spacecraft/iss-14.avif"], "texAlt": [null, null, null, null, null, null, null, null, null, null, null, null, null, null, null]}, "tiangong": {"mats": [{"name": "offwhite", "col": [0.965, 0.961, 0.933], "spec": 0.35, "tex": -1}, {"name": "gray", "col": [0.471, 0.49, 0.51], "spec": 0.237, "tex": -1}, {"name": "blue2", "col": [0.098, 0.149, 0.282], "spec": 0.9, "tex": -1}, {"name": "gold", "col": [0.588, 0.471, 0.275], "spec": 0.7, "tex": -1}, {"name": "dark", "col": [0.118, 0.133, 0.149], "spec": 0.237, "tex": -1}, {"name": "black", "col": [0.059, 0.059, 0.071], "spec": 0.237, "tex": -1}, {"name": "red", "col": [0.745, 0.118, 0.118], "spec": 0.237, "tex": -1}], "credit": null, "pack": "assets/models/satellites/tiangong.bin", "tex": [], "texAlt": []}, "colibri": {"mats": [{"name": "Material.001", "col": [0.711, 0.711, 0.711], "spec": 0.593, "tex": 0}], "nozzles": [[-2.24, -0.79, -11.2, 0.44], [2.24, -0.79, -11.2, 0.44]], "credit": {"title": "space craft", "author": "CyberDriger83", "url": "https://sketchfab.com/3d-models/space-craft-f7026b90bf9b44c99c15f7afc87bcdd3"}, "pack": "assets/models/ships/colibri.bin", "tex": ["assets/textures/spacecraft/colibri-0.avif"], "texAlt": [null]}, "ranger": {"mats": [{"name": "material", "col": [0.973, 0.973, 0.973], "spec": 0.15, "tex": 0}, {"name": ".002", "col": [1, 1, 1], "spec": 0.429, "tex": 1}, {"name": ".003", "col": [1, 1, 1], "spec": 0.15, "tex": 2}, {"name": ".006", "col": [1, 1, 1], "spec": 0.15, "tex": 3}, {"name": ".007", "col": [1, 1, 1], "spec": 0.15, "tex": 4}, {"name": ".008", "col": [1, 1, 1], "spec": 0.15, "tex": 5}, {"name": ".010", "col": [0.05, 0.07, 0.09], "spec": 1.2, "tex": -1}], "nozzles": [[-0.95, 0.29, -12.2, 0.62], [-0.02, -1.24, -12.2, 0.62], [1.03, 0.21, -12.2, 0.62]], "credit": {"title": "Interstellar Ranger One", "author": "Max Vizell", "url": "https://sketchfab.com/3d-models/interstellar-ranger-one-77c63df2062d4fd9863cc64711450c6f"}, "pack": "assets/models/ships/ranger.bin", "tex": ["assets/textures/spacecraft/ranger-0.jpg", "assets/textures/spacecraft/ranger-1.avif", "assets/textures/spacecraft/ranger-2.avif", "assets/textures/spacecraft/ranger-3.avif", "assets/textures/spacecraft/ranger-4.avif", "assets/textures/spacecraft/ranger-5.avif"], "texAlt": [null, null, null, null, null, null]}, "atlas": {"mats": [{"name": "Engine", "col": [1, 1, 1], "spec": 0.15, "tex": 0, "etex": 1, "estr": 0.8}, {"name": "Head", "col": [1, 1, 1], "spec": 0.15, "tex": 2, "etex": 3, "estr": 0.8}, {"name": "Middle", "col": [1, 1, 1], "spec": 0.15, "tex": 4, "etex": 3, "estr": 0.8}, {"name": "Back", "col": [1, 1, 1], "spec": 0.15, "tex": 5}, {"name": "Weapon", "col": [1, 1, 1], "spec": 0.15, "tex": 6}], "nozzles": [[-10.95, -0.03, -18.7, 1.17], [-6.04, 0.62, -18.7, 1.17], [6.04, 0.62, -18.7, 1.17], [10.98, -0.07, -18.7, 1.17]], "credit": {"title": "COBR-A", "author": "Alastoir De Faruh", "url": "https://sketchfab.com/3d-models/cobr-a-14588b866bb5461c853947e4b55da750"}, "pack": "assets/models/ships/atlas.bin", "tex": ["assets/textures/spacecraft/atlas-0.avif", "assets/textures/spacecraft/atlas-1.jpg", "assets/textures/spacecraft/atlas-2.avif", "assets/textures/spacecraft/atlas-3.jpg", "assets/textures/spacecraft/atlas-4.avif", "assets/textures/spacecraft/atlas-5.avif", "assets/textures/spacecraft/atlas-6.avif"], "texAlt": [null, null, null, null, null, null, null]}, "odisea": {"mats": [{"name": "Stuff", "col": [1, 1, 1], "spec": 0.15, "tex": 0, "etex": 1, "estr": 2.27}, {"name": "Main", "col": [1, 1, 1], "spec": 0.15, "tex": 2, "etex": 3, "estr": 2.5}], "nozzles": [[-3.49, 0.26, -70.2, 2.8], [-0.06, 0.35, -70.2, 2.8], [3.51, 0.24, -70.2, 2.8]], "credit": {"title": "SpaceShip AV-007", "author": "Muffin_one", "url": "https://sketchfab.com/3d-models/spaceship-av-007-1cd9c6349e714d5ead082177dc9e1250"}, "pack": "assets/models/ships/odisea.bin", "tex": ["assets/textures/spacecraft/odisea-0.avif", "assets/textures/spacecraft/odisea-1.jpg", "assets/textures/spacecraft/odisea-2.avif", "assets/textures/spacecraft/odisea-3.jpg"], "texAlt": [null, null, null, null]}, "parker": {"mats": [{"name": "foil_silver", "col": [0.654, 0.654, 0.654], "spec": 0.74, "tex": -1}, {"name": "base_metal", "col": [0.531, 0.529, 0.508], "spec": 0.675, "tex": -1}, {"name": "black_matte", "col": [0.304, 0.304, 0.304], "spec": 0.393, "tex": -1}, {"name": "material", "col": [0.728, 0.728, 0.728], "spec": 0.12, "tex": -1}, {"name": "foil_antenna", "col": [0.654, 0.654, 0.654], "spec": 0.762, "tex": -1}, {"name": "material_0", "col": [0.005, 0.005, 0.005], "spec": 0.138, "tex": -1}], "tex": [], "texAlt": [], "credit": {"title": "Parker Solar Probe", "author": "salvey", "url": "https://sketchfab.com/3d-models/parker-solar-probe-665bf4ab9de8400e9486bcb337f98e2c", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "adaptado: orientación (escudo térmico hacia el Sol), escala, simplificación y formato"}, "pack": "assets/models/satellites/parker.bin"}, "newhorizons": {"mats": [{"name": "lamina_metalica", "col": [1, 1, 1], "spec": 0.95, "tex": 0}, {"name": "pintura_blanca", "col": [1, 1, 1], "spec": 0.3, "tex": 0}, {"name": "superficie_oscura", "col": [1, 1, 1], "spec": 0.15, "tex": 0}], "tex": ["assets/textures/spacecraft/newhorizons-0.webp"], "texAlt": [null], "credit": {"title": "New Horizons", "author": null, "url": null, "license": null, "licenseUrl": null, "note": "modelo adquirido por el creador del proyecto; información de atribución pendiente. Adaptado: orientación (antena hacia la Tierra), escala real, simplificación, color compuesto a partir de los mapas difuso y especular, y formato"}, "pack": "assets/models/satellites/newhorizons.bin"}, "rocinante": {"mats": [{"name": "00_erged", "col": [1, 1, 1], "spec": 0.15, "tex": 0}], "tex": ["assets/textures/spacecraft/rocinante-0.webp"], "texAlt": [null], "nozzles": [[-2.45, 2.14, -23.2, 1.26], [-1.99, -3.03, -23.2, 1.26], [1.99, -3.03, -23.2, 1.26], [2.45, 2.14, -23.2, 1.26]], "credit": {"title": "Rosinante", "author": "mohamedhussien", "url": "https://sketchfab.com/3d-models/rosinante-9fb2988e9ede4721bff2145a3d07b7fd", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "adaptado: orientación (proa hacia +Z), escala a 46 m, simplificación y formato"}, "pack": "assets/models/ships/rocinante.bin"}, "orville": {"mats": [{"name": "Material", "col": [1, 1, 1], "spec": 0.15, "tex": 0, "etex": 1, "estr": 0.8}], "tex": ["assets/textures/spacecraft/orville-0.webp", "assets/textures/spacecraft/orville-1.webp"], "texAlt": [null, null], "nozzles": [[-9.45, 5.89, -150.2, 6.0], [1.19, 4.21, -150.2, 6.0], [11.75, 5.71, -150.2, 6.0]], "credit": {"title": "USS Orville", "author": "HaughtyGrayAlien", "url": "https://sketchfab.com/3d-models/uss-orville-60170fd125b348349cedeab68f3cca32", "license": "CC BY 4.0", "licenseUrl": "http://creativecommons.org/licenses/by/4.0/", "note": "adaptado: orientación (proa hacia +Z), escala a 300 m y formato"}, "pack": "assets/models/ships/orville.bin"}};
['colibri', 'ranger', 'atlas', 'odisea'].forEach(id => { const c = MODEL_PACKS[id] && MODEL_PACKS[id].credit; if (c) MODEL_CREDITS[id] = { title: c.title, author: c.author, url: c.url, license: 'CC BY 4.0', licenseUrl: 'https://creativecommons.org/licenses/by/4.0/', note: 'adaptado: orientación, escala, texturas y formato' }; });
const PackedModels = {
async gunzip(b64) {
const bin = Uint8Array.from(atob(b64), ch => ch.charCodeAt(0));
return new Response(new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
},
decode(buf) {
const dv = new DataView(buf);
if (String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3)) !== 'SLR2') return null;
const nv = dv.getUint32(4, true), ni = dv.getUint32(8, true), step = dv.getFloat32(12, true), ng = dv.getUint32(16, true);
let off = 20; const groups = [];
for (let g = 0; g < ng; g++) { groups.push({ first: dv.getUint32(off, true), count: dv.getUint32(off + 4, true), mat: dv.getUint32(off + 8, true) }); off += 12; }
const al = () => { off += (4 - (off % 4)) % 4; };
al(); const q = new Int16Array(buf, off, nv * 3); off += nv * 6; al();
const nm = new Int8Array(buf, off, nv * 4); off += nv * 4;
const uh = new Uint16Array(buf, off, nv * 2); off += nv * 4; al();
const dl = new Int32Array(buf, off, ni);
return { nv, ni, step, groups, q, nm, uh, dl };
},
h2f(h) {
const s = (h & 0x8000) ? -1 : 1, e = (h >> 10) & 0x1f, f = h & 0x3ff;
if (e === 0) return s * Math.pow(2, -14) * (f / 1024);
if (e === 31) return f ? NaN : s * Infinity;
return s * Math.pow(2, e - 15) * (1 + f / 1024);
},
async loadPack(name) {
const M = MODEL_PACKS[name];
if (!M || typeof DecompressionStream === 'undefined') return false;
const D = this.decode(await this.gunzip(M.pack)); if (!D) return false;
const { nv, ni, step, groups, q, nm, uh, dl } = D;
const vmat = new Uint8Array(nv);
const idx = new Uint32Array(ni); let a = 0; for (let i = 0; i < ni; i++) { a += dl[i]; idx[i] = a; }
groups.forEach(g => { for (let i = g.first; i < g.first + g.count; i++) vmat[idx[i]] = g.mat; });
const pos = new Float32Array(nv * 3), nrm = new Float32Array(nv * 3), col = new Float32Array(nv * 4), uv = new Float32Array(nv * 2);
let R = 0;
for (let i = 0; i < nv; i++) {
const x = q[i * 3] * step, y = q[i * 3 + 1] * step, z = q[i * 3 + 2] * step;
pos[i * 3] = x; pos[i * 3 + 1] = y; pos[i * 3 + 2] = z; R = Math.max(R, Math.hypot(x, y, z));
const nx = nm[i * 4], ny = nm[i * 4 + 1], nz = nm[i * 4 + 2], l = Math.hypot(nx, ny, nz) || 1;
nrm[i * 3] = nx / l; nrm[i * 3 + 1] = ny / l; nrm[i * 3 + 2] = nz / l;
const m = M.mats[vmat[i]] || { col: [0.6, 0.6, 0.6], spec: 0.3 };
col[i * 4] = m.col[0]; col[i * 4 + 1] = m.col[1]; col[i * 4 + 2] = m.col[2]; col[i * 4 + 3] = m.spec;
uv[i * 2] = this.h2f(uh[i * 2]); uv[i * 2 + 1] = this.h2f(uh[i * 2 + 1]);
}
for (let i = 0; i < nv * 3; i++) pos[i] /= R;
const srcs = await Promise.all(M.tex.map((u, i) => Assets.image(u, M.texAlt && M.texAlt[i])));
if (srcs.some(s => !s)) return false;      // sin un formato compatible: se conserva el modelo simplificado
const texs = await Promise.all(srcs.map(u => GLX.texFromURI(u).catch(() => null)));
const mesh = GLX.mesh({ a_pos: { data: pos, size: 3 }, a_nrm: { data: nrm, size: 3 }, a_col: { data: col, size: 4 }, a_uv: { data: uv, size: 2 } }, idx);
mesh.R = R; mesh.tris = ni / 3;
mesh.groups = groups.map(g => { const mt = M.mats[g.mat] || {}; return { first: g.first, count: g.count, tex: mt.tex >= 0 ? texs[mt.tex] : null, etex: mt.etex >= 0 ? texs[mt.etex] : null, estr: mt.estr || 0 }; });
CraftModels.setMesh(name, mesh);
return true;
},
async loadAll() {
for (const name of Object.keys(MODEL_PACKS)) {
try {
if (!(await this.loadPack(name))) continue;
const mesh = CraftModels.mesh(name);
World.rb.filter(rb => rb.isCraft && rb.def.model === name).forEach(rb => {
rb.lodMesh = CraftModels.mesh(name + '-lite'); rb.mesh = mesh; rb.R = rb.def.R = mesh.R / 1000;
});
World.hd[name] = true;
if (MODEL_PACKS[name].nozzles) { SHIP_NOZZLES[name] = MODEL_PACKS[name].nozzles; if (Flight.state === 'hangar') Flight.renderHangar(); }
if (UI.sel && UI.sel.def.model === name) UI.renderInfo(UI.sel);
} catch (e) { console.warn('Modelo detallado no disponible: ' + name, e); }
await new Promise(r => setTimeout(r, 30));
}
},
async loadJWST() {
if (!JWST_PACK || JWST_PACK.length < 20 || typeof DecompressionStream === 'undefined') return false;
const bin = await Assets.bytes(JWST_PACK);
const buf = await new Response(new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
const dv = new DataView(buf);
if (String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3)) !== 'SLRM') return false;
const nv = dv.getUint32(4, true), ni = dv.getUint32(8, true), step = dv.getFloat32(12, true);
let off = 28;
const q = new Int16Array(buf, off, nv * 3); off += nv * 6;
const nm = new Int8Array(buf, off, nv * 4); off += nv * 4; off += (4 - (off % 4)) % 4;
const dl = new Int32Array(buf, off, ni);
const pos = new Float32Array(nv * 3), nrm = new Float32Array(nv * 3), col = new Float32Array(nv * 4), idx = new Uint32Array(ni);
let R = 0;
for (let i = 0; i < nv; i++) {
const x = q[i * 3] * step, y = q[i * 3 + 1] * step, z = q[i * 3 + 2] * step;
pos[i * 3] = x; pos[i * 3 + 1] = y; pos[i * 3 + 2] = z;
R = Math.max(R, Math.hypot(x, y, z));
const nx = nm[i * 4], ny = nm[i * 4 + 1], nz = nm[i * 4 + 2], l = Math.hypot(nx, ny, nz) || 1;
nrm[i * 3] = nx / l; nrm[i * 3 + 1] = ny / l; nrm[i * 3 + 2] = nz / l;
const c = JWST_PALETTE[nm[i * 4 + 3] & 255] || [0.6, 0.6, 0.6, 0.3];
col[i * 4] = c[0]; col[i * 4 + 1] = c[1]; col[i * 4 + 2] = c[2]; col[i * 4 + 3] = c[3];
}
for (let i = 0; i < nv * 3; i++) pos[i] /= R;
let acc = 0; for (let i = 0; i < ni; i++) { acc += dl[i]; idx[i] = acc; }
const mesh = GLX.mesh({ a_pos: { data: pos, size: 3 }, a_nrm: { data: nrm, size: 3 }, a_col: { data: col, size: 4 } }, idx);
mesh.R = R; mesh.tris = ni / 3;
CraftModels.setMesh('jwst', mesh);
BODY.jwst.R = R / 1000;
return true;
},
};
const ASSET_MANIFEST = {"music": {"webm": "assets/audio/metamorphosis.webm", "mp3": "assets/audio/metamorphosis.mp3"}, "bodyMeshes": {"fobos": {"url": "assets/models/moons/fobos.bin", "k": 1.26255}, "deimos": {"url": "assets/models/moons/deimos.bin", "k": 1.37497}, "haumea": {"url": "assets/models/dwarfs/haumea.bin", "k": 1.44255}}};
const SOLARIS_BUILD = '2026.10.06-0855-web';
console.info('SOLARIS · versión ' + SOLARIS_BUILD);
const Assets = {
cache: new Map(), busyN: new Map(),
isURL: s => typeof s === 'string' && /^(\.\/)?assets\//.test(s),
async bytes(s) {
if (this.isURL(s)) {
const r = await fetch(s); if (!r.ok) throw new Error(s + ' → ' + r.status);
return new Uint8Array(await r.arrayBuffer());
}
return Uint8Array.from(atob(s), ch => ch.charCodeAt(0));
},
async gunzip(s) {
const bin = await this.bytes(s);
return new Response(new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
},
async image(primary, alt) {
if (!primary) return alt || null;
const avif = /^data:image\/avif|\.avif$/.test(primary);
if (!avif || await GLX.avifOK()) return primary;
return alt || null;
},
once(key, fn) { if (!this.cache.has(key)) this.cache.set(key, fn().catch(e => { console.warn('Recurso no disponible: ' + key, e); this.cache.delete(key + ':ok'); return false; })); return this.cache.get(key); },
busy(key, label) {
this.busyN.set(key, label); this.paintBusy();
},
done(key) { this.busyN.delete(key); this.paintBusy(); },
paintBusy() {
const el = document.getElementById('asset-load'); if (!el) return;
const labels = [...this.busyN.values()];
el.classList.toggle('on', labels.length > 0);
if (labels.length) el.querySelector('span').textContent = labels[labels.length - 1];
},
};
Object.assign(Assets, {
visual(rb) {
const id = rb.id, m = rb.def.model, st = s => s === 'ok' ? 'loaded' : s === 'fail' ? 'failed' : s === 'loading' ? 'loading' : 'available';
if (rb.isCraft && m && (MODEL_PACKS[m] || (m === 'jwst' && JWST_PACK && JWST_PACK.length > 20)))
return { type: 'mesh', status: st(PackedModels.state[m]), source: m };
if (ASSET_MANIFEST.bodyMeshes && ASSET_MANIFEST.bodyMeshes[id])
return { type: 'mesh', status: st(PackedModels.bodyState[id]), source: ASSET_MANIFEST.bodyMeshes[id].url };
if (typeof PLANET_TEX !== 'undefined' && PLANET_TEX[id]) {
const s = World.texState && World.texState[id];
return { type: s === 'fail' ? 'standard' : 'texture', status: st(s), source: PLANET_TEX[id].map };
}
return { type: 'standard', status: 'loaded', source: null };
},
changed() { clearTimeout(this.vbT); this.vbT = setTimeout(() => UI.refreshBadges && UI.refreshBadges(), 60); },
});
const PREV_CREDITS = [];
for (const id in PLANET_TEX) {
const c = PLANET_TEX[id].credit; if (!c || !BODY[id]) continue;
if (PLANET_TEX[id].prevCredit) PREV_CREDITS.push({ id, credit: PLANET_TEX[id].prevCredit });
if (PLANET_TEX[id].rings && BODY[id] && !BODY[id].vis.rings) BODY[id].vis.rings = PLANET_TEX[id].rings;   // anillo incluido en el modelo
if (BODY[id].credit && BODY[id].credit.url !== c.url) PREV_CREDITS.push({ id, credit: BODY[id].credit });
BODY[id].credit = c;
}
for (const name in MODEL_PACKS) if (MODEL_PACKS[name].credit && !MODEL_CREDITS[name]) { MODEL_CREDITS[name] = MODEL_PACKS[name].credit; if (BODY[name]) BODY[name].credit = MODEL_CREDITS[name]; }
Object.assign(PackedModels, {
state: {}, bodyState: {},
gunzip(s) { return Assets.gunzip(s); },
request(name, visible) {
if (this.state[name] === 'ok' || this.state[name] === 'fail') return Promise.resolve(this.state[name] === 'ok');
if (visible) Assets.busy('m:' + name, 'Cargando modelo…');
return Assets.once('model:' + name, async () => {
this.state[name] = 'loading'; Assets.changed();
let ok = false;
try { ok = name === 'jwst' ? await this.loadJWST() : await this.loadPack(name); } catch (e) { console.warn('Modelo detallado no disponible: ' + name, e); }
this.state[name] = ok ? 'ok' : 'fail'; Assets.changed();
if (ok) this.apply(name);
Assets.done('m:' + name);
return ok;
}).finally(() => Assets.done('m:' + name));
},
apply(name) {
const mesh = CraftModels.mesh(name);
World.rb.filter(rb => rb.isCraft && rb.def.model === name).forEach(rb => {
rb.mesh = mesh; rb.lodMesh = CraftModels.mesh(name + '-lite');
if (mesh.R) rb.R = rb.def.R = mesh.R / 1000;
});
World.hd[name] = true; if (name === 'jwst') World.jwstHD = true;
if (MODEL_PACKS[name] && MODEL_PACKS[name].nozzles) { SHIP_NOZZLES[name] = MODEL_PACKS[name].nozzles; }
if (Flight.ship && Flight.ship.id === name && Flight.syncModel) Flight.syncModel();
if (Flight.state === 'hangar') Flight.renderHangar();
if (UI.sel && UI.sel.def.model === name) UI.renderInfo(UI.sel);
},
loadBody(id) {
const E = ASSET_MANIFEST.bodyMeshes && ASSET_MANIFEST.bodyMeshes[id];
if (!E || typeof DecompressionStream === 'undefined') return Promise.resolve(false);
this.bodyState[id] = this.bodyState[id] || 'loading'; Assets.changed();
return Assets.once('body:' + id, async () => {
let D = null; try { D = this.decode(await Assets.gunzip(E.url)); } catch (e) { console.warn('Forma no disponible: ' + id, e); }
if (!D) { this.bodyState[id] = 'fail'; Assets.changed(); return false; }
const { nv, ni, step, q, nm, dl } = D;
const pos = new Float32Array(nv * 3), nrm = new Float32Array(nv * 3), dir = new Float32Array(nv * 3); let R = 0;
for (let i = 0; i < nv * 3; i++) pos[i] = q[i] * step;
for (let i = 0; i < nv; i++) R = Math.max(R, Math.hypot(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]));
for (let i = 0; i < nv; i++) {
const x = pos[i * 3] / R * E.k, y = pos[i * 3 + 1] / R * E.k, z = pos[i * 3 + 2] / R * E.k, l = Math.hypot(x, y, z) || 1;
pos[i * 3] = x; pos[i * 3 + 1] = y; pos[i * 3 + 2] = z;                 // radio medio = 1 (escala visual de SOLARIS)
dir[i * 3] = x / l; dir[i * 3 + 1] = y / l; dir[i * 3 + 2] = z / l;
const nx = nm[i * 4], ny = nm[i * 4 + 1], nz = nm[i * 4 + 2], n = Math.hypot(nx, ny, nz) || 1;
nrm[i * 3] = nx / n; nrm[i * 3 + 1] = ny / n; nrm[i * 3 + 2] = nz / n;
}
const idx = new Uint32Array(ni); let a = 0; for (let i = 0; i < ni; i++) { a += dl[i]; idx[i] = a; }
const rb = World.byId[id]; if (!rb) return false;
const old = rb.mesh;
rb.mesh = GLX.mesh({ a_pos: { data: pos, size: 3 }, a_nrm: { data: nrm, size: 3 }, a_dir: { data: dir, size: 3 } }, idx);
rb.irregular = true; rb.extent = E.k; World.hd[id] = true; this.bodyState[id] = 'ok'; Assets.changed();   // extensión real: radio máximo / radio medio
if (UI.sel === rb) UI.renderInfo(rb);
if (old && old !== rb.mesh && GLX.freeMesh) GLX.freeMesh(old);           // libera la forma procedural anterior
return true;
});
},
async preloadShips(first) {
const order = [first, ...SHIPS.map(s => s.id).filter(id => id !== first)].filter(id => MODEL_PACKS[id]);
for (const id of order) { await this.request(id, id === first); await new Promise(r => setTimeout(r, 60)); }
},
tick() {
for (const rb of World.rb) {
if (!rb.isCraft || !rb.def.model) continue;
const name = rb.def.model; if (this.state[name]) continue;
if (!(MODEL_PACKS[name] || name === 'jwst')) continue;
const focus = UI.sel === rb || UI.tourRB === rb || (Flight.on && Flight.target === rb);
if ((rb.drawOn && rb.proj.rpx > 9) || focus) this.request(name, focus);
}
},
});
const STARS = {
betelgeuse: ['Betelgeuse', 5.919, 7.407, 0.5], rigel: ['Rigel', 5.242, -8.202, 0.13], bellatrix: ['Bellatrix', 5.419, 6.35, 1.64],
saiph: ['Saiph', 5.796, -9.67, 2.07], alnitak: ['Alnitak', 5.679, -1.943, 1.77], alnilam: ['Alnilam', 5.604, -1.202, 1.69],
mintaka: ['Mintaka', 5.533, -0.299, 2.23], meissa: ['Meissa', 5.586, 9.934, 3.33],
dubhe: ['Dubhe', 11.062, 61.751, 1.79], merak: ['Merak', 11.031, 56.382, 2.37], phecda: ['Phecda', 11.897, 53.695, 2.44],
megrez: ['Megrez', 12.257, 57.033, 3.31], alioth: ['Alioth', 12.9, 55.96, 1.77], mizar: ['Mizar', 13.399, 54.925, 2.23], alkaid: ['Alkaid', 13.792, 49.313, 1.86],
caph: ['Caph', 0.153, 59.15, 2.28], schedar: ['Schedar', 0.675, 56.537, 2.24], navi: ['Navi', 0.945, 60.717, 2.4], ruchbah: ['Ruchbah', 1.43, 60.235, 2.68], segin: ['Segin', 1.907, 63.67, 3.37],
deneb: ['Deneb', 20.69, 45.28, 1.25], sadr: ['Sadr', 20.37, 40.257, 2.23], albireo: ['Albireo', 19.512, 27.96, 3.05], gienah: ['Gienah', 20.77, 33.97, 2.48], deltacyg: ['δ Cyg', 19.75, 45.131, 2.87],
vega: ['Vega', 18.616, 38.784, 0.03], sheliak: ['Sheliak', 18.835, 33.363, 3.5], sulafat: ['Sulafat', 18.982, 32.69, 3.25], deltalyr: ['δ Lyr', 18.908, 36.899, 4.3], zetalyr: ['ζ Lyr', 18.746, 37.605, 4.3],
altair: ['Altair', 19.846, 8.868, 0.76], tarazed: ['Tarazed', 19.771, 10.613, 2.72], alshain: ['Alshain', 19.922, 6.407, 3.71], deltaaql: ['δ Aql', 19.425, 3.115, 3.36], lambdaaql: ['λ Aql', 19.104, -4.882, 3.44], zetaaql: ['ζ Aql', 19.09, 13.863, 2.99], thetaaql: ['θ Aql', 20.188, -0.821, 3.26],
antares: ['Antares', 16.49, -26.432, 0.96], graffias: ['Graffias', 16.091, -19.806, 2.62], dschubba: ['Dschubba', 16.006, -22.622, 2.29], pisco: ['π Sco', 15.981, -26.114, 2.89], sigmasco: ['σ Sco', 16.353, -25.593, 2.88],
tausco: ['τ Sco', 16.598, -28.216, 2.82], epssco: ['ε Sco', 16.836, -34.293, 2.29], musco: ['μ Sco', 16.864, -38.047, 3.0], zetasco: ['ζ Sco', 16.91, -42.362, 3.62], etasco: ['η Sco', 17.203, -43.239, 3.33],
sargas: ['Sargas', 17.622, -42.998, 1.86], iotasco: ['ι Sco', 17.793, -40.127, 2.99], kappasco: ['κ Sco', 17.708, -39.03, 2.39], shaula: ['Shaula', 17.56, -37.104, 1.62], lesath: ['Lesath', 17.513, -37.296, 2.7],
kausaus: ['Kaus Australis', 18.403, -34.385, 1.85], nunki: ['Nunki', 18.921, -26.297, 2.05], ascella: ['Ascella', 19.044, -29.88, 2.6], kausmed: ['Kaus Media', 18.35, -29.828, 2.7],
kausbor: ['Kaus Borealis', 18.466, -25.422, 2.81], alnasl: ['Alnasl', 18.097, -30.424, 2.98], phisgr: ['φ Sgr', 18.761, -26.991, 3.17], tausgr: ['τ Sgr', 19.116, -27.67, 3.32],
regulus: ['Régulo', 10.14, 11.967, 1.35], denebola: ['Denébola', 11.818, 14.572, 2.14], algieba: ['Algieba', 10.333, 19.842, 2.08], zosma: ['Zosma', 11.235, 20.524, 2.56], chertan: ['Chertan', 11.237, 15.43, 3.33],
etaleo: ['η Leo', 10.122, 16.763, 3.48], adhafera: ['Adhafera', 10.278, 23.417, 3.44], muleo: ['μ Leo', 9.879, 26.007, 3.88], epsleo: ['ε Leo', 9.764, 23.774, 2.98],
castor: ['Cástor', 7.577, 31.888, 1.58], pollux: ['Pólux', 7.755, 28.026, 1.14], alhena: ['Alhena', 6.629, 16.399, 1.93], mebsuta: ['Mebsuta', 6.732, 25.131, 3.0], tejat: ['Tejat', 6.383, 22.514, 2.87], wasat: ['Wasat', 7.335, 21.982, 3.53], alzirr: ['Alzirr', 6.755, 12.896, 3.35],
aldebaran: ['Aldebarán', 4.599, 16.509, 0.86], elnath: ['Elnath', 5.438, 28.608, 1.65], zetatau: ['ζ Tau', 5.627, 21.143, 3.0], gammatau: ['γ Tau', 4.33, 15.628, 3.65], deltatau: ['δ Tau', 4.382, 17.543, 3.76],
epstau: ['ε Tau', 4.477, 19.18, 3.53], thetatau: ['θ Tau', 4.478, 15.871, 3.4], lambdatau: ['λ Tau', 4.011, 12.49, 3.4], alcyone: ['Alcíone (Pléyades)', 3.791, 24.105, 2.87],
sirius: ['Sirio', 6.752, -16.716, -1.46], mirzam: ['Mirzam', 6.378, -17.956, 1.98], adhara: ['Adhara', 6.977, -28.972, 1.5], wezen: ['Wezen', 7.14, -26.393, 1.84], aludra: ['Aludra', 7.402, -29.303, 2.45], muliphein: ['Muliphein', 7.063, -15.633, 4.1],
acrux: ['Ácrux', 12.443, -63.099, 0.77], mimosa: ['Mimosa', 12.795, -59.689, 1.25], gacrux: ['Gácrux', 12.519, -57.113, 1.64], deltacru: ['δ Cru', 12.253, -58.749, 2.79],
rigilkent: ['Rigil Kentaurus', 14.66, -60.834, -0.27], hadar: ['Hadar', 14.064, -60.373, 0.61],
arcturus: ['Arturo', 14.261, 19.182, -0.05], izar: ['Izar', 14.75, 27.074, 2.37], muphrid: ['Muphrid', 13.911, 18.398, 2.68], seginus: ['Seginus', 14.535, 38.308, 3.03], deltaboo: ['δ Boo', 15.258, 33.315, 3.47], nekkar: ['Nekkar', 15.032, 40.39, 3.5], rhoboo: ['ρ Boo', 14.531, 30.371, 3.58],
spica: ['Espiga', 13.42, -11.161, 0.97], porrima: ['Porrima', 12.694, -1.449, 2.74], vindemiatrix: ['Vindemiatrix', 13.036, 10.959, 2.83], deltavir: ['δ Vir', 12.927, 3.397, 3.38], betavir: ['β Vir', 11.845, 1.765, 3.6], zetavir: ['ζ Vir', 13.578, -0.596, 3.37],
polaris: ['Polaris', 2.53, 89.264, 1.98], kochab: ['Kochab', 14.845, 74.156, 2.08], pherkad: ['Pherkad', 15.346, 71.834, 3.0], yildun: ['Yildun', 17.537, 86.586, 4.36], epsumi: ['ε UMi', 16.766, 82.037, 4.2], zetaumi: ['ζ UMi', 15.734, 77.795, 4.3], etaumi: ['η UMi', 16.292, 75.755, 4.95],
markab: ['Markab', 23.079, 15.205, 2.49], scheat: ['Scheat', 23.063, 28.083, 2.42], algenib: ['Algenib', 0.221, 15.184, 2.83], alpheratz: ['Alpheratz', 0.14, 29.091, 2.06], enif: ['Enif', 21.736, 9.875, 2.38], homam: ['Homam', 22.691, 10.831, 3.4],
mirach: ['Mirach', 1.162, 35.621, 2.05], almach: ['Almach', 2.065, 42.33, 2.1], deltaand: ['δ And', 0.655, 30.861, 3.27],
mirfak: ['Mirfak', 3.405, 49.861, 1.79], algol: ['Algol', 3.136, 40.956, 2.1], gammaper: ['γ Per', 3.08, 53.506, 2.93], deltaper: ['δ Per', 3.715, 47.788, 3.01], epsper: ['ε Per', 3.964, 40.01, 2.89], zetaper: ['ζ Per', 3.902, 31.884, 2.85],
capella: ['Capella', 5.278, 45.998, 0.08], menkalinan: ['Menkalinan', 5.992, 44.948, 1.9], thetaaur: ['θ Aur', 5.995, 37.213, 2.62], iotaaur: ['ι Aur', 4.95, 33.166, 2.69], epsaur: ['ε Aur', 5.033, 43.823, 3.0],
procyon: ['Proción', 7.655, 5.225, 0.34], gomeisa: ['Gomeisa', 7.453, 8.289, 2.89],
canopus: ['Canopo', 6.399, -52.696, -0.74], achernar: ['Achernar', 1.629, -57.237, 0.46], fomalhaut: ['Fomalhaut', 22.961, -29.622, 1.16],
};
const CONSTELLATIONS = [
{ id: 'orion', name: 'Orión', lines: [['meissa', 'betelgeuse'], ['meissa', 'bellatrix'], ['betelgeuse', 'alnitak'], ['bellatrix', 'mintaka'], ['alnitak', 'alnilam'], ['alnilam', 'mintaka'], ['alnitak', 'saiph'], ['mintaka', 'rigel']] },
{ id: 'umayor', name: 'Osa Mayor', lines: [['dubhe', 'merak'], ['merak', 'phecda'], ['phecda', 'megrez'], ['megrez', 'dubhe'], ['megrez', 'alioth'], ['alioth', 'mizar'], ['mizar', 'alkaid']] },
{ id: 'umenor', name: 'Osa Menor', lines: [['polaris', 'yildun'], ['yildun', 'epsumi'], ['epsumi', 'zetaumi'], ['zetaumi', 'kochab'], ['kochab', 'pherkad'], ['pherkad', 'etaumi'], ['etaumi', 'zetaumi']] },
{ id: 'casiopea', name: 'Casiopea', lines: [['caph', 'schedar'], ['schedar', 'navi'], ['navi', 'ruchbah'], ['ruchbah', 'segin']] },
{ id: 'cisne', name: 'Cisne', lines: [['deneb', 'sadr'], ['sadr', 'albireo'], ['deltacyg', 'sadr'], ['sadr', 'gienah']] },
{ id: 'lira', name: 'Lira', lines: [['vega', 'zetalyr'], ['zetalyr', 'deltalyr'], ['deltalyr', 'sulafat'], ['sulafat', 'sheliak'], ['sheliak', 'zetalyr']] },
{ id: 'aguila', name: 'Águila', lines: [['tarazed', 'altair'], ['altair', 'alshain'], ['altair', 'deltaaql'], ['deltaaql', 'lambdaaql'], ['zetaaql', 'deltaaql'], ['alshain', 'thetaaql']] },
{ id: 'escorpion', name: 'Escorpión', lines: [['graffias', 'dschubba'], ['dschubba', 'pisco'], ['dschubba', 'sigmasco'], ['sigmasco', 'antares'], ['antares', 'tausco'], ['tausco', 'epssco'], ['epssco', 'musco'], ['musco', 'zetasco'], ['zetasco', 'etasco'], ['etasco', 'sargas'], ['sargas', 'iotasco'], ['iotasco', 'kappasco'], ['kappasco', 'shaula'], ['shaula', 'lesath']] },
{ id: 'sagitario', name: 'Sagitario', lines: [['alnasl', 'kausmed'], ['kausmed', 'kausaus'], ['kausaus', 'alnasl'], ['kausmed', 'kausbor'], ['kausbor', 'phisgr'], ['phisgr', 'kausmed'], ['phisgr', 'nunki'], ['nunki', 'tausgr'], ['tausgr', 'ascella'], ['ascella', 'phisgr'], ['ascella', 'kausaus']] },
{ id: 'leo', name: 'Leo', lines: [['regulus', 'etaleo'], ['etaleo', 'algieba'], ['algieba', 'adhafera'], ['adhafera', 'muleo'], ['muleo', 'epsleo'], ['algieba', 'zosma'], ['zosma', 'denebola'], ['denebola', 'chertan'], ['chertan', 'regulus'], ['zosma', 'chertan']] },
{ id: 'geminis', name: 'Géminis', lines: [['castor', 'pollux'], ['castor', 'mebsuta'], ['mebsuta', 'tejat'], ['pollux', 'wasat'], ['wasat', 'alhena'], ['wasat', 'alzirr']] },
{ id: 'tauro', name: 'Tauro', lines: [['elnath', 'epstau'], ['epstau', 'deltatau'], ['deltatau', 'gammatau'], ['gammatau', 'thetatau'], ['thetatau', 'aldebaran'], ['aldebaran', 'zetatau'], ['gammatau', 'lambdatau']] },
{ id: 'cmayor', name: 'Can Mayor', lines: [['mirzam', 'sirius'], ['sirius', 'muliphein'], ['sirius', 'wezen'], ['wezen', 'adhara'], ['wezen', 'aludra']] },
{ id: 'cmenor', name: 'Can Menor', lines: [['procyon', 'gomeisa']] },
{ id: 'cruz', name: 'Cruz del Sur', lines: [['acrux', 'gacrux'], ['mimosa', 'deltacru']] },
{ id: 'centauro', name: 'Centauro (punteros)', lines: [['rigilkent', 'hadar']] },
{ id: 'boyero', name: 'Boyero', lines: [['arcturus', 'izar'], ['izar', 'deltaboo'], ['deltaboo', 'nekkar'], ['nekkar', 'seginus'], ['seginus', 'rhoboo'], ['rhoboo', 'arcturus'], ['arcturus', 'muphrid']] },
{ id: 'virgo', name: 'Virgo', lines: [['betavir', 'porrima'], ['porrima', 'deltavir'], ['deltavir', 'vindemiatrix'], ['porrima', 'spica'], ['deltavir', 'zetavir'], ['zetavir', 'spica']] },
{ id: 'pegaso', name: 'Pegaso y Andrómeda', lines: [['markab', 'scheat'], ['scheat', 'alpheratz'], ['alpheratz', 'algenib'], ['algenib', 'markab'], ['markab', 'homam'], ['homam', 'enif'], ['alpheratz', 'deltaand'], ['deltaand', 'mirach'], ['mirach', 'almach']] },
{ id: 'perseo', name: 'Perseo', lines: [['gammaper', 'mirfak'], ['mirfak', 'deltaper'], ['deltaper', 'epsper'], ['epsper', 'zetaper'], ['mirfak', 'algol']] },
{ id: 'auriga', name: 'Auriga', lines: [['capella', 'menkalinan'], ['menkalinan', 'thetaaur'], ['thetaaur', 'elnath'], ['elnath', 'iotaaur'], ['iotaaur', 'epsaur'], ['epsaur', 'capella']] },
];
const DEEP_SKY = [
['Nebulosa de Orión', 5.588, -5.39, 1.4, [1.0, 0.45, 0.55], 0.9, 1, 0],
['Nebulosa de la Laguna', 18.06, -24.38, 1.1, [1.0, 0.5, 0.55], 0.7, 1, 0],
['Nebulosa de Carina', 10.75, -59.87, 2.2, [1.0, 0.45, 0.5], 0.8, 1, 0],
['Nebulosa Roseta', 6.53, 4.95, 1.2, [1.0, 0.4, 0.45], 0.45, 1, 0],
['Nebulosa Norteamérica', 20.98, 44.33, 2.2, [1.0, 0.42, 0.45], 0.5, 1, 0],
['Galaxia de Andrómeda', 0.712, 41.27, 2.2, [1.0, 0.92, 0.82], 0.9, 3.2, 38],
['Gran Nube de Magallanes', 5.39, -69.76, 5.0, [0.86, 0.9, 1.0], 0.85, 1.3, 170],
['Pequeña Nube de Magallanes', 0.877, -72.8, 2.6, [0.88, 0.92, 1.0], 0.65, 1.6, 45],
['Galaxia del Triángulo', 1.564, 30.66, 0.9, [0.9, 0.92, 1.0], 0.35, 1.6, 23],
['Pléyades', 3.791, 24.1, 1.0, [0.6, 0.75, 1.0], 0.55, 1, 0],
['Saco de Carbón', 12.83, -62.5, 3.0, [0, 0, 0], -1.0, 1.3, 30],
];
const EARTH_LAND = [
'-168,65.6 -166,68.3 -163,70 -156.5,71.3 -152,70.8 -145,70.1 -140,69.6 -134,69.6 -128,70.1 -120,68.8 -114,68.2 -108,68.6 -101,67.8 -96,68.2 -94,71.8 -90,68.6 -85,69.8 -82,68 -81.5,66.5 -86,66 -88,64.5 -93,63 -94.5,61 -94.7,59 -93,58 -90,57.2 -85.5,55.3 -82.2,55.1 -81,52.5 -79.5,51.3 -78.9,52.5 -78.7,55 -77,56.3 -77.2,58.3 -78.3,60 -77.8,62.4 -74,62.2 -72,61.8 -69.5,59.2 -67.5,58.5 -65,60.3 -62,57.5 -60,55.5 -57,53 -56,51.6 -60,50.2 -64.5,49.2 -66.5,50.1 -70,47.5 -64.5,48.8 -65,47.5 -61,45.6 -63.5,44.6 -66,43.7 -66.5,45 -70,43.8 -70.6,42.6 -70,41.7 -74,40.5 -74.9,38.9 -76,37 -75.5,35.2 -78,33.9 -81,31.5 -80.5,28.5 -80,26 -80.4,25.2 -81.8,26.2 -82.6,27.8 -84,30 -86,30.4 -89,30.2 -89.2,29.1 -91,29.5 -94,29.6 -97,27.8 -97.5,25.5 -97.7,22.3 -96,19.2 -94.5,18.2 -91.5,18.6 -90.3,21.1 -87.5,21.5 -86.8,21 -87.5,18.5 -88.3,16.5 -86,15.9 -83.3,15.1 -83.7,11 -83,10 -81,8.9 -79.6,9.6 -77.4,8.7 -77.9,7.2 -78.5,8.3 -80,7.3 -81.7,8.1 -83.6,8.4 -85.7,10.3 -87.5,12.9 -90,13.8 -92.2,14.5 -94.5,16.2 -96.5,15.7 -99.9,16.8 -102,17.9 -105.6,20.4 -105.3,21.6 -106.4,23.2 -108,25.5 -109.9,27.6 -111,29.5 -114.7,31.7 -114.2,29.5 -112.5,27 -110.3,24.2 -109.4,23 -112,24.8 -112.1,25.9 -114.1,27.7 -115,28 -116,30 -117.1,32.5 -118.5,34 -120.6,34.6 -122.5,37.5 -124.2,40.3 -124.4,43 -124,46.2 -124.7,48.4 -127,50.5 -128.2,52.5 -130.5,54.7 -133,57 -136,58.3 -140,59.7 -144,60.1 -148,60.5 -152,59.2 -156,57.5 -162,55 -164.5,54.5 -161,56.2 -158,58.6 -162,58.6 -164.8,60.4 -165.5,61.8 -164.5,63.1 -161,64.4 -166.5,64.6',
'-77.4,8.7 -76,10 -75.5,10.8 -74.2,11.3 -72.2,12.2 -71.6,10.7 -70,12 -68,10.6 -64.5,10.3 -62,10.7 -61,9 -60,8.4 -58.3,6.8 -55,6 -52,5 -51,4.2 -50,1.8 -49.5,0 -48,-1 -44.5,-2.3 -41,-2.9 -38.5,-3.7 -35.2,-5.4 -34.8,-7.5 -35.3,-9.6 -37.2,-11.4 -38.6,-13 -39,-17.7 -40.3,-20.5 -41.6,-22.6 -43.2,-23 -44.8,-23.4 -48.5,-26 -48.6,-28.3 -50.5,-30.8 -52.5,-33.5 -53.8,-34.6 -56,-34.9 -57.5,-35.5 -57.2,-36.3 -56.8,-37.2 -57.6,-38.2 -62,-39 -62.3,-40.6 -65,-41 -64,-42.3 -65,-44.8 -67.5,-46.3 -66,-47.8 -68.3,-50.1 -69,-51.6 -68.4,-52.4 -68.6,-54.8 -65.3,-54.9 -67.3,-55.9 -70,-55.2 -72.5,-53.5 -74.8,-52 -75.5,-48.5 -74.5,-46 -74,-43.5 -73.8,-41.5 -73.6,-39 -73.2,-37 -72,-34 -71.6,-32.5 -71.4,-29.5 -70.5,-25.5 -70.2,-20 -70.3,-18.3 -72.5,-16.8 -75.2,-15.3 -76.3,-13.6 -77.2,-12 -78.8,-8.8 -80.8,-6 -81.3,-4.5 -80.3,-3.3 -80.9,-2.2 -80.4,-0.4 -80,1 -78.9,1.6 -77.4,4.1 -77.3,6.5 -77.9,7.2',
'-5.9,35.8 -2,35.1 1,36.5 3,36.8 6.5,37.1 9.8,37.3 11.1,36.9 10.3,35.2 10.9,34.1 10.1,33.5 11.5,33.1 15.2,32.3 17,31 19.5,30.3 20,32.1 22,32.9 25,31.8 29,30.9 31,31.6 32.3,31.3 32.5,29.9 33.5,27.5 35.5,23.9 37.2,21 38.5,18 39.7,15.5 41.5,13.9 43.3,12.6 43.4,11.7 44.5,10.4 48,11.2 51.2,11.8 51,10.4 50.8,8.3 49,6 47.5,4.5 46,2 43.5,-0.5 42,-1.2 40.5,-2.5 39.5,-4.7 39.2,-6.5 39.5,-8 40.4,-10.5 40.6,-14.5 40.5,-15.5 37,-17.8 35,-20 35.5,-22 35.6,-23.7 33,-25.8 32.9,-26.3 32.4,-28.6 31,-29.9 29,-32.1 27,-33.7 25.6,-34 22,-34.2 20,-34.8 18.4,-34.3 17.9,-32.5 16.5,-28.6 15.1,-26.6 14.5,-22.9 13.3,-20.5 11.8,-17.2 12.3,-13.5 13.4,-10.5 13.1,-8.8 12.3,-6 11.8,-4.5 9.6,-2.5 9.3,-0.7 9.5,2 9.7,3.6 8.5,4.5 7,4.4 5.7,4.3 4,6.4 1.5,6.2 -2.1,4.75 -4,5.2 -7.5,4.4 -9.5,5.5 -11.5,6.9 -13.2,8.5 -14.5,10.7 -16.7,12.3 -17.2,14.7 -16.5,16.2 -16.1,19.6 -17.05,20.8 -15,24 -14.5,26.1 -12.9,27.9 -11,28.8 -9.8,30.3 -9.3,32.5 -6.8,34.1',
'49.3,-12 50.5,-15.5 49.5,-17.5 48.5,-21 47,-25.2 45,-25.5 43.6,-22 44.4,-17.5 44,-16.5 46.3,-15.7 48,-14',
'-5.6,36 -6.4,36.8 -7.4,37.2 -8.9,37 -9,38.7 -8.7,41.2 -9.3,43 -8,43.7 -5,43.5 -1.8,43.4 -1.3,44.5 -1.2,46.2 -2.3,47.2 -4.7,48 -4.8,48.4 -3,48.8 -1.6,48.7 -1.6,49.7 0.2,49.5 1.6,50.2 1.9,51 3.5,51.4 4.5,52.5 4.7,53 6.8,53.5 8.5,53.6 8.6,54.9 8.1,55.6 8.6,57.1 10.6,57.7 10.3,56.2 10,55 10.9,54 13.5,54.3 16,54.5 18.6,54.6 21,55 21,56.8 24,57.2 23.5,59 26,59.5 28,59.7 30,59.9 29,60.2 25,60.2 22,60.3 21.4,61.5 21.5,63 25.3,65 24.5,65.8 22,65.5 21,64.5 17.8,62.5 17.1,61 18.7,60 17,58.7 16.5,57 14.2,55.4 12.8,55.6 12.5,56.5 11,58.9 10.6,59.9 8,58.1 5.6,58.9 5,60.5 5,62 8,63.5 11,64.9 13,66.5 15.5,68.5 18.5,69.8 23.5,70.9 25.8,71.1 28.5,70.9 31,70.3 33,69.4 36,69.1 41,67.6 40.5,66 44,66.5 44,68.5 46,68.3 53.5,68.3 58,68.7 60,69.8 66,69 68,72.5 70,73 72.5,72.5 73.5,68.5 76,72 80,73.5 86,74.5 95,76 104,77.7 110,76.5 113.5,73.6 126,73.4 130,71 139,72 146,72.3 152,70.8 159.7,70.8 166,69.5 170,70 176,69.8 180,68.9 180,65 177,64.5 178.5,62.5 174,61.8 170,60 164,59.8 163,57.8 162.3,56 160,53 156.7,51 156,52.5 155.6,55 156.8,57.8 158,58 160,60.6 155,59.3 151,59.5 143,59.3 140.5,58.5 137,54.5 140.5,53.5 141.4,52.2 141,49 140.5,48 138.5,46.5 135.1,43.5 131.9,43 130.7,42.3 129.5,41 129.4,37 129.3,35.2 127.5,34.6 126.3,34.6 126.5,36.5 126.2,37.7 125,39.6 124.3,40 121.5,39 121.2,40.8 119.5,39.9 117.7,39 118.5,38 120.5,37.8 122.6,37.4 120.3,36 119.2,35 120.8,32.2 121.9,30.9 122,29.9 121.5,28.5 120,26.5 118.5,24.5 116.5,23 114.2,22.3 112,21.6 110.5,21 110.2,20.3 108.5,21.6 106.7,20.5 105.8,19 107,16.8 108.8,15.3 109.4,12.5 109,11.4 107,10.4 105,8.6 104.8,10.4 103,11.5 102,12.5 100.9,13.3 100,12 99.2,10 100.3,8.3 101.2,6.9 102.5,6 103.4,4 103.5,2 104.2,1.4 103.3,1.4 101.3,2.8 100.4,4.2 100.2,6 98.3,8 98.6,10 97.7,15.4 97.4,16.5 94.8,16 94.2,18.5 93,20 92.2,21.5 91.8,22.3 90.5,22 88.8,21.7 86.9,21 86,19.8 84.5,18.5 82,16.5 80.3,15.5 80.2,13 79.8,10.3 77.5,8.1 76.3,9.8 75,12.5 73.5,16 72.8,19 72.6,21.3 70.5,20.8 69,22.4 68.5,23.5 67,24.8 64,25.3 61.5,25.1 58,25.6 57.2,26.5 56.3,27.2 54,26.8 51.5,27.9 50.3,29.9 48.5,30 48,29.5 48.8,27.8 50.2,26.2 51.6,25 51.6,26.1 51.2,24.5 54.5,24.2 56.3,26.2 56.4,24.9 58.5,23.6 59.8,22.5 58.5,20.5 57.5,19 55,17 52.2,15.7 49,14.2 45,12.8 43.5,12.7 42.7,15.5 41.5,17.5 39.2,21.5 38,24 36.5,25.9 35,28 35,29.5 34.6,28.1 34.25,27.73 32.6,29.9 32.3,31.3 34.2,31.3 34.9,32.8 35.5,34 35.9,35.5 36.2,36.6 34.5,36.8 32.5,36.1 30.5,36.3 28.5,36.7 27.3,37.8 26.5,39.5 26.2,40.1 27,40.9 26,40.8 24,40.5 22.9,40.6 23.5,39 24,38 23,36.5 21.7,36.8 21.1,38.3 20.2,39.6 19.4,41.8 18.5,42.5 16,43.5 14,45 13.6,45.6 12.3,45.4 12.6,44.1 13.6,43.5 16,41.9 17.2,40.9 18.5,40.1 17,39.4 16,38.1 15.6,40 14.3,40.8 12.6,41.5 11.1,42.4 10.3,43.6 8.9,44.4 7,43.7 5,43.3 3.2,43.1 3.2,42 0.9,41 -0.3,39.5 0.2,38.7 -0.7,37.6 -2.1,36.7 -4.4,36.7',
'-180,68.9 -176,67.5 -172,66.2 -169.7,66 -171.5,64.5 -173,64.3 -176,65.2 -180,65',
'114.1,-21.8 116.7,-20.6 118.8,-20.3 121,-19.5 122.2,-18 122.2,-17 123.5,-16.3 125,-14.6 126.9,-13.9 128.2,-15 129.5,-14.9 130.2,-13 130.8,-12.4 132.6,-11.5 135,-12.2 136.8,-12.2 136.5,-13.3 135.5,-15 137.7,-16.3 140.5,-17.6 141.6,-15 141.5,-12.5 142.5,-10.7 143.5,-12.5 144.5,-14.3 145.4,-16 146,-18.5 148.8,-20.3 150.8,-22.5 153.1,-25.2 153.6,-28.6 153,-31 151.3,-33.8 150,-37.3 148,-37.8 146.2,-39.1 144.5,-38.2 143.5,-38.8 141,-38.1 139.7,-37.2 138.5,-35.6 137.8,-33.2 136,-34.8 135.2,-33.5 133,-32 131,-31.5 128,-32.2 124,-33 121.5,-33.8 118,-35 115,-34.3 115.7,-32 115.1,-30 114,-27 113.5,-24.5',
'144.7,-40.7 148.3,-40.9 148,-43.2 146.5,-43.6 145.2,-42.2',
'172.7,-34.4 174.5,-35.5 175.9,-37.4 178.5,-37.7 177.9,-39.2 176.8,-40 175.2,-41.6 174.6,-41.3 174.7,-39.3 173.8,-39.2 174.5,-37.3 173,-35',
'172.7,-40.5 174.2,-41.7 173.3,-43 172.7,-43.8 171.3,-44.4 170.8,-45.9 169,-46.6 166.5,-46 166.6,-45.2 168.3,-44 170.5,-42.9 172,-41.4',
'131,-1.4 134,-0.9 135.5,-3.3 137.9,-1.5 141,-2.6 144.5,-3.8 146,-5.5 147.9,-6.2 150.8,-10.2 149.5,-10.4 147.5,-10.1 146,-8.1 144,-7.7 142.7,-9.3 141,-9.1 138.8,-8.4 137.5,-5.1 135,-4.4 133,-4 132,-2.8',
'109,1.5 109.6,2 111.5,2.5 113,3.2 115.4,5 117,7 117.5,6.5 119.2,5.2 118.1,4.3 118,2.5 119,1 117.8,0.5 117.5,-1.3 116.5,-2.5 116,-3.9 114.6,-4.1 113,-3.2 111.6,-3.4 110.2,-2.9 110,-1.3 109.1,0',
'95.3,5.6 97.5,5.2 100.4,2.3 103.7,-0.9 104.7,-2.5 106,-3.3 105.8,-5.9 104.5,-5.9 102.3,-4 100.3,-1 98.7,1.7 96.4,3.9',
'105.2,-6.8 106.1,-5.9 108.3,-6.3 110.4,-6.9 112.6,-6.9 114.5,-7.7 114.4,-8.7 111,-8.2 108,-7.8 106.4,-7.4',
'119.4,-5.5 119.6,-3.5 118.8,-2.6 119.8,0.2 120.8,1.3 124.3,1.5 125.2,1.6 124,0.4 121,0.5 120.2,-1 121.6,-1.9 123.3,-1 122.8,-4.5 121.5,-4.8 120.4,-2.9 120.4,-5.6',
'120.6,18.5 122.2,18.5 122,16 121.6,14.2 124,13.2 123.2,13.8 120.6,14.4 120,16',
'122,7 123.5,8.6 125.4,9.7 126.5,7.3 125.4,5.6 124,6.3',
'130.9,34 132.5,35.4 135.4,35.6 136.8,37.3 139.2,38 140,40 141.4,41.4 142,39.5 141,37 140.8,35.7 139.8,34.9 138.8,34.6 137,34.6 135.8,33.4 135.1,34.3 133,34.4 131,33.9',
'140.2,41.5 141.2,41.8 143.3,42 145.5,43.3 144.5,44 141.9,45.5 141.4,43.2 140.4,43.3 139.9,42.5',
'129.8,33.4 131.2,33.9 132,33 131.6,31.6 130.6,31 130.2,31.6 129.6,32.8',
'132.5,33.2 133.9,34.4 134.7,34.1 134.3,33.4 133,32.8',
'121.5,25.3 121.9,24.6 120.9,22 120.1,23 120.1,24',
'80.2,9.8 81.9,7.4 81.6,6.4 80.5,5.9 79.8,7',
'108.6,19.2 110.1,20.1 111,19.6 110.5,18.5 109.5,18.2 108.7,18.5',
'142,46 143.5,46.3 143,49 144.7,49 143,52.5 142.8,54.3 142.2,54.2 141.7,52.5 142.1,50 141.9,47.5',
'-5.7,50.1 -3.5,50.4 0,50.8 1.4,51.2 1.7,52.7 0.3,53.5 -0.1,54.5 -1.6,55.6 -2.1,57.1 -1.8,57.6 -3.3,58.6 -5,58.6 -6.2,56.7 -5.6,55.4 -4.8,54.9 -3,54.1 -3,53.4 -4.7,53.3 -4.2,52.3 -5.2,51.7 -4.2,51.2 -4.6,50.7',
'-6.2,53.3 -6,54.4 -7.3,55.3 -8.5,54.6 -10.2,54.2 -9.9,53.4 -10.3,51.9 -9.5,51.5 -8,51.8 -6.4,52.2',
'-22.7,64.1 -24,65.5 -22.5,66.4 -16,66.5 -14.5,66.1 -13.6,65.1 -15,64.3 -18,63.4 -21,63.8',
'12.4,38 15.6,38.3 15.1,36.7 12.6,37.6', '8.4,41.2 9.8,41 9.6,39.2 8.4,38.9', '8.6,42.9 9.5,43 9.4,41.4 8.7,41.8',
'23.5,35.6 26.3,35.3 25.8,34.9 23.6,35.2', '32.3,35 34.6,35.7 33.8,34.8 32.6,34.7',
'-73,78.5 -66,81 -55,82.2 -40,83.5 -25,83 -18,81.5 -12,81.5 -18,79 -19.5,76 -21,73.5 -22,70.5 -25,69.5 -32,68.2 -37.5,65.6 -40.5,64 -43,60 -44,59.8 -48,61 -50.5,63.8 -53,66.5 -54,69.5 -51.5,70.5 -55,71.5 -56.5,74.5 -60,76 -66,76 -71.5,77.5',
'11,78.8 16,80 22,80.3 27,79.5 21,78 17,76.6 14,77.5',
'52,71.5 56,74 62,76.2 68.5,76.9 64,75.4 57,72.5 55,70.8',
'-63,66.5 -61.8,66.8 -68.5,70.5 -73,72 -80,73.7 -89.5,72.8 -85,70.5 -82,69.8 -78.5,67 -74,64.5 -71,62.8 -65.5,62 -64.5,63.5',
'-118,69 -102,68.5 -101,70 -105,73 -114,73.5 -119,71.5',
'-80,76.5 -73,78.5 -64,82 -75,83 -90,81.5 -95,78.5 -89,76.3',
'-80,74.5 -91,74.6 -92,76.3 -82,76.6',
'-125,71.9 -120.5,71.4 -116,73.5 -119,74.6 -124.5,74.3',
'-59.4,47.6 -55.5,51.6 -53,49.5 -52.7,47.5 -53.5,46.7 -56,47.5',
'-84.9,21.9 -82,23.2 -79.5,22.8 -77,21.3 -74.1,20.2 -77.5,19.9 -80,21.7 -82.5,21.6',
'-74.4,18.4 -72.8,19.9 -70,19.8 -68.3,18.6 -70,18.2 -72,18.2',
'-180,-78 -165,-78.2 -155,-77.5 -150,-76.5 -140,-75 -130,-74.3 -120,-73.8 -110,-74 -100,-73 -90,-72.8 -80,-73.3 -75,-71 -68,-67 -65,-65 -60,-63.5 -57,-63.3 -59,-64.5 -61,-67 -63,-69 -61,-72 -60,-75 -55,-77.5 -45,-78 -36,-78 -30,-77 -25,-75 -20,-73.5 -10,-71 0,-70 10,-70 20,-69.8 30,-69.5 40,-68.8 50,-66.5 60,-67.3 70,-69.5 75,-69 80,-67.5 90,-66.5 100,-65.7 110,-66 120,-66.5 130,-66 140,-66.5 150,-68.5 160,-70 165,-71.5 168,-74 165,-77.5 170,-78 180,-78 180,-90 -180,-90',
];
const EARTH_WATER = [
'-94.8,59 -93,58 -90,57.2 -85.5,55.3 -82.2,55.1 -81,52.5 -79.5,51.3 -78.9,52.5 -78.7,55 -77,56.3 -77.2,58.3 -78.3,60 -77.5,62.3 -80.5,63.1 -84,63.5 -87,63.7 -90.5,63.5 -92.5,62 -94.3,60.8',
'28,41.2 28.6,43.4 29.9,45.3 30.7,46.5 32.5,45.4 33.5,44.4 35,44.8 35.5,46.3 39.3,47.2 38,45.3 37.8,44.7 39.7,43.5 41.6,41.6 39.7,41 37,41.1 35,42 33,41.9 31.5,41.2 29.1,41.2',
'47,44.9 49,46.5 51.5,47 53,46.6 53,45.2 51.3,44.5 52.8,41.8 53,40 53.9,38.9 53.9,37.3 51.5,36.8 49.5,37.5 48.8,38.7 49.5,40.4 48.6,41.8 47.5,43 46.8,44.5',
'-92,46.8 -88,48.5 -84.5,46.5 -84.8,45.8 -87.5,44 -87.8,41.7 -86.5,42.5 -85,46 -88,46.5',
'-84,45.9 -80.5,45 -79.7,43.3 -76.5,43.6 -79,42.8 -82.4,41.7 -83.5,42 -82.5,43 -82.4,45.3',
];
const EARTH_RANGES = [
'-75,10 -77,3 -78,-2 -76,-10 -70,-16 -68,-22 -69.5,-30 -70.5,-36 -71.5,-42 -73,-49 -72,-53',
'-152,63 -140,61 -130,58 -125,54 -118,50 -112,45 -108,40 -106,36 -105,32 -106,28 -103,22 -98,19 -92,16',
'70,36 75,35.5 80,32 85,28.5 90,28 95,29 98,31', '78,35 84,34 90,33.5 96,33', '6,45 10,46.5 14,47', '40,43.5 48,41',
'45,35 52,30 57,27', '-9,31 0,34 9,36', '60,50 59,58 60,65 65,68', '38,8 37,0 35,-8', '146,-16 148,-24 150,-30 148,-36',
'7,61 14,66 20,69', '-84,34 -78,39 -72,44', '-121,49 -121,40 -118,36', '100,28 103,24', '44,42 46,40',
];
const SETTINGS_KEY = 'solaris.settings', SETTINGS_VERSION = 1;
const GFX_PRESETS = {
low: { scale: 0.75, aa: false, shadows: 1, textures: 'low', distance: 'low', lod: 'low', particles: 'low', bloom: false, sunGlow: true, atmospheres: true, nebulae: false, dust: false, sky: 'basic', fps: 60 },
medium: { scale: 1, aa: true, shadows: 2, textures: 'medium', distance: 'medium', lod: 'medium', particles: 'medium', bloom: true, sunGlow: true, atmospheres: true, nebulae: true, dust: true, sky: 'high', fps: 60 },
high: { scale: 1, aa: true, shadows: 3, textures: 'high', distance: 'high', lod: 'high', particles: 'high', bloom: true, sunGlow: true, atmospheres: true, nebulae: true, dust: true, sky: 'high', fps: 0 },
ultra: { scale: 1.25, aa: true, shadows: 3, textures: 'ultra', distance: 'max', lod: 'high', particles: 'high', bloom: true, sunGlow: true, atmospheres: true, nebulae: true, dust: true, sky: 'max', fps: 0 },
};
const GFX_KEYS = Object.keys(GFX_PRESETS.high);
const KEY_ACTIONS = {
explore: [
['search', '/', 'Buscar'], ['center', 'c', 'Centrar la cámara'], ['follow', 'f', 'Seguir al objeto'], ['explore', 'e', 'Explorar la selección'],
['overview', 'h', 'Vista general'], ['sun', 's', 'Regresar al Sol'], ['zoomIn', '+', 'Acercar'], ['zoomOut', '-', 'Alejar'],
['rotL', 'arrowleft', 'Girar a la izquierda'], ['rotR', 'arrowright', 'Girar a la derecha'], ['rotU', 'arrowup', 'Girar hacia arriba'], ['rotD', 'arrowdown', 'Girar hacia abajo'],
['pause', 'space', 'Pausar o reanudar el tiempo'], ['faster', '.', 'Tiempo más rápido'], ['slower', ',', 'Tiempo más lento'], ['now', 'r', 'Fecha y hora actuales'],
['names', 'l', 'Mostrar u ocultar nombres'], ['orbits', 'o', 'Mostrar u ocultar órbitas'], ['constel', 'k', 'Constelaciones'], ['minimap', 'm', 'Minimapa'],
['cine', 'v', 'Modo cine'], ['tours', 't', 'Recorridos guiados'], ['pilot', 'p', 'Pilotar una nave'],
],
flight: [
['thrUp', 'w', 'Acelerar (subir empuje)'], ['thrDown', 's', 'Reducir empuje'], ['thr0', 'x', 'Motores en reposo'], ['brake', 'space', 'Frenar'],
['left', 'a', 'RCS: desplazar a la izquierda'], ['right', 'd', 'RCS: desplazar a la derecha'], ['up', 'r', 'RCS: subir'], ['down', 'f', 'RCS: bajar'],
['pitchUp', 'arrowup', 'Cabeceo hacia arriba'], ['pitchDown', 'arrowdown', 'Cabeceo hacia abajo'], ['yawL', 'arrowleft', 'Guiñada a la izquierda'], ['yawR', 'arrowright', 'Guiñada a la derecha'],
['rollL', 'q', 'Alabeo a la izquierda'], ['rollR', 'e', 'Alabeo a la derecha'], ['cam', 'c', 'Cambiar cámara'], ['target', 't', 'Seleccionar objetivo'],
['match', 'b', 'Igualar velocidad'], ['aim', 'g', 'Orientar hacia el objetivo'], ['intercept', 'i', 'Piloto automático'], ['hold', 'h', 'Mantener distancia'],
['orbit', 'o', 'Entrar o salir de órbita'], ['autobrake', 'n', 'Freno de aproximación'], ['stab', 'z', 'Estabilización'], ['mode', 'v', 'Modo asistido o simulación'],
['warpUp', '.', 'Más aceleración temporal'], ['warpDown', ',', 'Menos aceleración temporal'], ['sound', 'm', 'Sonido de cabina'],
],
};
const KEY_RESERVED = new Set(['escape', 'home', 'tab', 'enter', '?', 'shift', 'control', 'alt', 'meta', 'altgraph', 'capslock', 'f1', 'f5', 'f11', 'f12', 'pageup', 'pagedown', 'contextmenu', 'dead', 'unidentified', ...'0123456789']);
const Settings = {
state: null, _saveT: 0, listeners: [],
defaults() {
const rm = !!(window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches);
const ctl = {}; for (const c in KEY_ACTIONS) { ctl[c] = {}; KEY_ACTIONS[c].forEach(([a, k]) => ctl[c][a] = k); }
return {
v: SETTINGS_VERSION,
general: { startup: 'intro', confirmExit: true, tips: true },
audio: { master: 0.9, music: 0.30, ui: 0.6, nav: 0.6, ship: 0.8, hover: true, muted: false, musicOn: true },
graphics: Object.assign({ auto: true, preset: 'high', ambient: 0.02 }, GFX_PRESETS.high),
interface: { uiScale: 1, panelOpacity: 0.6, info: 'advanced', hud: 'full', anim: 'full', minimap: true },
controls: ctl,
simulation: { timeSpeed: 'real', startDate: 'now', lastJD: null, scale: 'visual', rotCap: true,
layers: { orbits: true, moonOrbits: true, moons: true, asteroids: true, kuiper: true, craft: true, constLines: false, constNames: false, namesPlanets: true, namesMoons: true, namesSmall: true, namesCraft: true, hz: false }, hzModel: 'conservador', sciPreset: null, plannerMode: 'simple' },
language: 'es',
accessibility: { reducedMotion: rm, contrast: false, text: 'normal', reduceFlashes: false },
firstRun: true,
};
},
merge(base, over) {
for (const k in over) {
if (!(k in base)) continue;                         // descarta claves desconocidas u obsoletas
if (base[k] && typeof base[k] === 'object' && !Array.isArray(base[k]) && over[k] && typeof over[k] === 'object') this.merge(base[k], over[k]);
else if (typeof over[k] === typeof base[k] || base[k] === null) base[k] = over[k];
}
return base;
},
load() {
const st = this.defaults();
let saved = null;
try { saved = JSON.parse(localStorage.getItem(SETTINGS_KEY) || 'null'); } catch (e) {}
if (saved && typeof saved === 'object') {
if ((saved.v || 0) <= SETTINGS_VERSION) this.merge(st, saved);    // versiones futuras: se ignoran con seguridad
st.firstRun = false;
} else {
try { const m = JSON.parse(localStorage.getItem('solaris.music') || 'null'); if (m) { st.audio.musicOn = m.on !== false; if (typeof m.vol === 'number') st.audio.music = m.vol; } } catch (e) {}
try { const a = JSON.parse(localStorage.getItem('solaris.audio') || 'null'); if (a) ['master', 'ui', 'ship', 'muted', 'hover'].forEach(k => { if (k in a) st.audio[k] = a[k]; }); if (a && 'ui' in a) st.audio.nav = a.ui; } catch (e) {}
}
st.v = SETTINGS_VERSION; st.language = 'es';
this.state = st;
return st;
},
save() {
clearTimeout(this._saveT); this._dirty = true;
this._saveT = setTimeout(() => this.flush(), 200);
},
flush() {
if (!this._dirty || !this.state) return;
clearTimeout(this._saveT); this._dirty = false;
try { localStorage.setItem(SETTINGS_KEY, JSON.stringify(Object.assign({}, this.state, { firstRun: false }))); } catch (e) {}
},
get(path) { return path.split('.').reduce((o, k) => o && o[k], this.state); },
set(path, value, opts) {
const ks = path.split('.'), last = ks.pop(), obj = ks.reduce((o, k) => o[k], this.state);
if (obj[last] === value && !(opts && opts.force)) return;
obj[last] = value;
if (ks[0] === 'graphics' && GFX_KEYS.includes(last) && !(opts && opts.fromPreset)) { this.state.graphics.preset = 'custom'; this.state.graphics.auto = false; }
this.save();
this.listeners.forEach(f => f(path, value));
},
on(f) { this.listeners.push(f); },
applyPreset(name) {
const P = GFX_PRESETS[name]; if (!P) return;
for (const k in P) this.set('graphics.' + k, P[k], { fromPreset: true });
this.state.graphics.preset = name; this.save();
this.listeners.forEach(f => f('graphics.preset', name));
},
reset() {
const keep = this.state.simulation.lastJD;
this.state = this.defaults(); this.state.firstRun = false; this.state.simulation.lastJD = keep;
if (this.state.graphics.auto) Object.assign(this.state.graphics, GFX_PRESETS[this.detectQuality()], { preset: this.detectQuality() });
this.save();
this.listeners.forEach(f => f('*', null));
},
resetControls() {
const d = this.defaults().controls; this.state.controls = d; this.save();
this.listeners.forEach(f => f('controls', d));
},
detectQuality(gl) {
gl = gl || (window.GLX && GLX.gl);
let r = '';
try { const ext = gl && gl.getExtension('WEBGL_debug_renderer_info'); r = (ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl ? gl.getParameter(gl.RENDERER) : '') || ''; } catch (e) {}
r = r.toLowerCase();
let level = 2;                                                         // 0 bajo · 1 medio · 2 alto
if (/swiftshader|llvmpipe|software|basic render|microsoft basic/.test(r)) level = 0;
else if (/intel|mali|adreno|powervr|videocore|uhd|iris/.test(r) && !/arc/.test(r)) level = 1;
const mem = navigator.deviceMemory, cores = navigator.hardwareConcurrency;
if (mem && mem <= 4) level = Math.min(level, 1);
if (cores && cores <= 4) level = Math.min(level, 1);
const px = screen.width * screen.height * Math.pow(Math.min(devicePixelRatio || 1, 2), 2);
if (px > 9e6 && level > 0) level -= 1;                                  // pantallas muy densas: un nivel menos
try { if (gl && gl.getParameter(gl.MAX_TEXTURE_SIZE) < 8192) level = Math.min(level, 1); } catch (e) {}
return ['low', 'medium', 'high'][level];
},
};
const Keys = {
norm(e) { const k = e.key; if (!k) return ''; if (k === ' ' || k === 'Spacebar') return 'space'; return k.length === 1 ? k.toLowerCase() : k.toLowerCase(); },
code(ctx, action) { return Settings.state.controls[ctx][action]; },
action(ctx, key) {
const m = Settings.state.controls[ctx];
for (const a in m) if (m[a] === key) return a;
const alias = ctx === 'explore' ? { '=': 'zoomIn', '_': 'zoomOut', ']': 'faster', '[': 'slower' } : { 'pageup': 'warpUp', 'pagedown': 'warpDown' };
const a = alias[key]; if (a && m[a] === KEY_ACTIONS[ctx].find(x => x[0] === a)[1]) return a;
return null;
},
label(k) {
const map = { space: 'Espacio', arrowup: '↑', arrowdown: '↓', arrowleft: '←', arrowright: '→', backspace: '⌫', delete: 'Supr', insert: 'Insert', end: 'Fin' };
return map[k] ? (window.I18N ? I18N.t(map[k]) : map[k]) : k.length === 1 ? k.toUpperCase() : k.charAt(0).toUpperCase() + k.slice(1);
},
describe(ctx, action) { const r = KEY_ACTIONS[ctx].find(x => x[0] === action); return r ? r[2] : action; },
assign(ctx, action, key) {
if (!key || KEY_RESERVED.has(key)) return { ok: false, reason: 'reserved' };
const m = Settings.state.controls[ctx], prev = m[action];
if (prev === key) return { ok: true };
const other = Object.keys(m).find(a => a !== action && m[a] === key);
m[action] = key; if (other) m[other] = prev;
Settings.save(); Settings.listeners.forEach(f => f('controls', null));
return { ok: true, swapped: other || null };
},
};
const Gfx = {
g() { return Settings.state.graphics; },
particles() { return { low: 0.35, medium: 0.65, high: 1 }[this.g().particles] || 1; },
lod() { return { low: [380, 80, 60], medium: [160, 30, 28], high: [90, 16, 12] }[this.g().lod] || [160, 30, 28]; },
occluders() { return [0, 1, 2, 4][this.g().shadows] == null ? 4 : [0, 1, 2, 4][this.g().shadows]; },
ringShadow() { return this.g().shadows >= 2; },
flash() { return Settings.state.accessibility.reduceFlashes ? 0.4 : 1; },
distPx() { return { low: 2.5, medium: 1.2, high: 0.6, max: 0.3 }[this.g().distance] || 0.6; },
};
const EVENT_CATS = { exploracion: 'Exploración espacial', astronomia: 'Fenómenos astronómicos' };
const EVENTS = [
{ date: '1957-10-04', cat: 'exploracion', title: 'Lanzamiento del Sputnik 1', desc: 'Primer satélite artificial puesto en órbita alrededor de la Tierra.', obj: 'tierra', mission: 'Sputnik 1' },
{ date: '1961-04-12', cat: 'exploracion', title: 'Primer vuelo espacial tripulado', desc: 'Yuri Gagarin completa una órbita alrededor de la Tierra a bordo de la Vostok 1.', obj: 'tierra', mission: 'Vostok 1' },
{ date: '1969-07-20', cat: 'exploracion', title: 'Primer alunizaje tripulado', desc: 'El módulo lunar del Apolo 11 se posa en el mar de la Tranquilidad.', obj: 'luna', mission: 'Apolo 11' },
{ date: '1976-07-20', cat: 'exploracion', title: 'Viking 1 aterriza en Marte', desc: 'Primer aterrizaje plenamente exitoso en la superficie marciana con envío prolongado de datos.', obj: 'marte', mission: 'Viking 1' },
{ date: '1977-08-20', cat: 'exploracion', title: 'Lanzamiento de la Voyager 2', desc: 'Inicio del «Gran Tour» por los planetas exteriores.', obj: 'voyager2', mission: 'Voyager 2' },
{ date: '1977-09-05', cat: 'exploracion', title: 'Lanzamiento de la Voyager 1', desc: 'La sonda que se convertiría en el objeto humano más lejano.', obj: 'voyager1', mission: 'Voyager 1' },
{ date: '1979-03-05', cat: 'exploracion', title: 'Voyager 1 sobrevuela Júpiter', desc: 'Descubre los volcanes activos de Ío.', obj: 'jupiter', mission: 'Voyager 1' },
{ date: '1979-07-09', cat: 'exploracion', title: 'Voyager 2 sobrevuela Júpiter', desc: 'Segundo sobrevuelo del sistema joviano por el programa Voyager.', obj: 'jupiter', mission: 'Voyager 2' },
{ date: '1980-11-12', cat: 'exploracion', title: 'Voyager 1 sobrevuela Saturno', desc: 'Incluye un sobrevuelo cercano de Titán.', obj: 'saturno', mission: 'Voyager 1' },
{ date: '1981-08-25', cat: 'exploracion', title: 'Voyager 2 sobrevuela Saturno', desc: 'La sonda continúa después hacia Urano.', obj: 'saturno', mission: 'Voyager 2' },
{ date: '1986-01-24', cat: 'exploracion', title: 'Voyager 2 sobrevuela Urano', desc: 'Única visita de una nave a Urano hasta la fecha.', obj: 'urano', mission: 'Voyager 2' },
{ date: '1989-08-25', cat: 'exploracion', title: 'Voyager 2 sobrevuela Neptuno', desc: 'Única visita de una nave a Neptuno; también estudia Tritón.', obj: 'neptuno', mission: 'Voyager 2' },
{ date: '1990-04-24', cat: 'exploracion', title: 'Lanzamiento del telescopio Hubble', desc: 'Telescopio espacial en órbita baja terrestre.', obj: 'hubble', mission: 'Hubble' },
{ date: '1995-12-02', cat: 'exploracion', title: 'Lanzamiento de SOHO', desc: 'Observatorio solar situado cerca del punto de Lagrange L1.', obj: 'soho', mission: 'SOHO' },
{ date: '1998-11-20', cat: 'exploracion', title: 'Primer módulo de la Estación Espacial Internacional', desc: 'Lanzamiento del módulo Zariá, inicio del ensamblaje de la EEI.', obj: 'iss', mission: 'EEI' },
{ date: '2004-07-01', cat: 'exploracion', title: 'Cassini llega a Saturno', desc: 'Inicio de trece años de estudio del sistema de Saturno.', obj: 'saturno', mission: 'Cassini-Huygens' },
{ date: '2005-01-14', cat: 'exploracion', title: 'Huygens aterriza en Titán', desc: 'Primer aterrizaje en un cuerpo del Sistema Solar exterior.', obj: 'titan', mission: 'Cassini-Huygens' },
{ date: '2006-01-19', cat: 'exploracion', title: 'Lanzamiento de New Horizons', desc: 'Rumbo a Plutón y al cinturón de Kuiper.', obj: 'newhorizons', mission: 'New Horizons' },
{ date: '2011-07-16', cat: 'exploracion', title: 'Dawn entra en órbita de Vesta', desc: 'Primera nave en orbitar un cuerpo del cinturón principal.', obj: 'vesta', mission: 'Dawn' },
{ date: '2012-08-06', cat: 'exploracion', title: 'Curiosity aterriza en Marte', desc: 'El vehículo explorador llega al cráter Gale.', obj: 'marte', mission: 'Mars Science Laboratory' },
{ date: '2012-08-25', cat: 'exploracion', title: 'Voyager 1 cruza la heliopausa', desc: 'Primera nave en entrar en el espacio interestelar.', obj: 'voyager1', mission: 'Voyager 1' },
{ date: '2015-03-06', cat: 'exploracion', title: 'Dawn entra en órbita de Ceres', desc: 'Primera nave en orbitar un planeta enano.', obj: 'ceres', mission: 'Dawn' },
{ date: '2015-07-14', cat: 'exploracion', title: 'New Horizons sobrevuela Plutón', desc: 'Primeras imágenes detalladas de Plutón y Caronte.', obj: 'pluton', mission: 'New Horizons' },
{ date: '2016-07-05', cat: 'exploracion', title: 'Juno entra en órbita de Júpiter', desc: 'Estudio del interior, la atmósfera y la magnetosfera de Júpiter.', obj: 'jupiter', mission: 'Juno' },
{ date: '2017-09-15', cat: 'exploracion', title: 'Fin de la misión Cassini', desc: 'La sonda se sumerge de forma controlada en la atmósfera de Saturno.', obj: 'saturno', mission: 'Cassini-Huygens' },
{ date: '2018-08-12', cat: 'exploracion', title: 'Lanzamiento de la sonda solar Parker', desc: 'Misión para estudiar la corona solar desde muy cerca.', obj: 'parker', mission: 'Parker Solar Probe' },
{ date: '2018-11-05', cat: 'exploracion', title: 'Voyager 2 cruza la heliopausa', desc: 'Segunda nave en alcanzar el espacio interestelar.', obj: 'voyager2', mission: 'Voyager 2' },
{ date: '2019-01-01', cat: 'exploracion', title: 'New Horizons sobrevuela Arrokoth', desc: 'El objeto más lejano visitado por una nave hasta entonces.', obj: 'arrokoth', mission: 'New Horizons' },
{ date: '2020-12-05', cat: 'exploracion', title: 'Hayabusa2 entrega muestras de Ryugu', desc: 'La cápsula con material del asteroide aterriza en Australia.', obj: 'ryugu', mission: 'Hayabusa2' },
{ date: '2021-02-18', cat: 'exploracion', title: 'Perseverance aterriza en Marte', desc: 'El vehículo explorador llega al cráter Jezero.', obj: 'marte', mission: 'Mars 2020' },
{ date: '2021-04-29', cat: 'exploracion', title: 'Lanzamiento del módulo Tianhe', desc: 'Módulo central de la estación espacial china Tiangong.', obj: 'tiangong', mission: 'Tiangong' },
{ date: '2021-12-25', cat: 'exploracion', title: 'Lanzamiento del telescopio James Webb', desc: 'Observatorio infrarrojo situado cerca del punto L2.', obj: 'jwst', mission: 'JWST' },
{ date: '2023-07-01', cat: 'exploracion', title: 'Lanzamiento de Euclid', desc: 'Telescopio para cartografiar la materia y la energía oscuras.', obj: 'euclid', mission: 'Euclid' },
{ date: '2023-09-24', cat: 'exploracion', title: 'OSIRIS-REx entrega muestras de Bennu', desc: 'La cápsula con material del asteroide aterriza en Utah.', obj: 'bennu', mission: 'OSIRIS-REx' },
{ date: '2024-12-24', cat: 'exploracion', title: 'Parker: máximo acercamiento al Sol', desc: 'La sonda pasa a unos 6.1 millones de km de la superficie solar.', obj: 'parker', mission: 'Parker Solar Probe' },
{ date: '1986-02-09', cat: 'astronomia', title: 'Perihelio del cometa Halley', desc: 'Máximo acercamiento al Sol en su último paso.', obj: 'halley' },
{ date: '1997-04-01', cat: 'astronomia', title: 'Perihelio del cometa Hale-Bopp', desc: 'Uno de los cometas más brillantes del siglo XX.', obj: 'halebopp' },
{ date: '2012-06-06', cat: 'astronomia', title: 'Tránsito de Venus', desc: 'Venus pasa por delante del disco solar visto desde la Tierra; el siguiente será en diciembre de 2117.', obj: 'venus' },
{ date: '2019-11-11', cat: 'astronomia', title: 'Tránsito de Mercurio', desc: 'Mercurio cruza el disco solar visto desde la Tierra.', obj: 'mercurio' },
{ date: '2024-04-08', cat: 'astronomia', title: 'Eclipse total de Sol', desc: 'Visible como total desde México, Estados Unidos y Canadá.', obj: 'luna' },
{ date: '2026-08-12', cat: 'astronomia', title: 'Eclipse total de Sol', desc: 'Franja de totalidad sobre Groenlandia, Islandia y España.', obj: 'luna' },
{ date: '2027-08-02', cat: 'astronomia', title: 'Eclipse total de Sol', desc: 'Franja de totalidad sobre el sur de España, el norte de África y Egipto.', obj: 'luna', pred: true },
{ date: '2029-04-13', cat: 'astronomia', title: 'Máximo acercamiento de Apofis', desc: 'El asteroide pasará a unos 32,000 km de la superficie terrestre, sin riesgo de impacto.', obj: 'apofis', pred: true },
{ date: '2032-11-13', cat: 'astronomia', title: 'Tránsito de Mercurio', desc: 'Próximo tránsito de Mercurio por delante del Sol.', obj: 'mercurio', pred: true },
{ date: '2061-07-28', cat: 'astronomia', title: 'Próximo perihelio del cometa Halley', desc: 'Fecha prevista de su regreso al Sistema Solar interior.', obj: 'halley', pred: true },
{ date: '2117-12-11', cat: 'astronomia', title: 'Tránsito de Venus', desc: 'Próximo tránsito de Venus por delante del Sol.', obj: 'venus', pred: true },
];
const HZ_MODELS = {
conservador: { name: 'Conservador', in: 0.99, out: 1.70, inName: 'pérdida de agua (efecto invernadero húmedo)', outName: 'máximo efecto invernadero',
ref: 'Kopparapu et al. (2013), ApJ 765, 131' },
optimista: { name: 'Optimista', in: 0.75, out: 1.77, calc: true, inName: '«Venus reciente»', outName: '«Marte primitivo»',
ref: 'Kopparapu et al. (2013): límites empíricos; distancias calculadas a partir del flujo efectivo (d = 1/√S, con S = 1.78 y 0.32)' },
clasico: { name: 'Clásico', in: 0.95, out: 1.67, inName: 'pérdida de agua', outName: 'máximo efecto invernadero', ref: 'Kasting et al. (1993)' },
};
const SYSTEMS = {
solar: { id: 'solar', name: 'Sistema Solar', sub: 'Nuestro sistema planetario', center: 'sol' },
gargantua: { id: 'gargantua', name: 'Sistema Gargantúa', sub: 'Sistema ficticio de Interstellar', center: 'gargantua', fiction: true,
H: 60,                                   // radio visual del horizonte de sucesos, en unidades de escena
disk: { in: 3.0, out: 11.0 },            // disco de acreción, en radios del horizonte (borde interior típico de un disco delgado: 3 radios)
maxDist: 60 * 240,
note: 'Sistema ficticio inspirado en Interstellar. Algunos conceptos representados están basados en fenómenos físicos reales.',
orbitsNote: 'Las órbitas del Sistema Gargantúa son una representación cinematográfica y educativa: la película no establece parámetros orbitales precisos.' },
};
const GARG_BODIES = [
{ id: 'gargantua', system: 'gargantua', name: 'Gargantúa', type: 'blackhole', R: 0, visR: 60, color: '#f0c27a', vis: { style: 'bh' }, extent: 3.4,
src: 'Película Interstellar (2014) y «The Science of Interstellar», de Kip Thorne (asesor científico de la película).',
info: {
desc: 'Agujero negro supermasivo en rápida rotación: el centro del sistema planetario que visita la misión Endurance.',
fiction: [['Tipo', 'Agujero negro supermasivo en rotación'], ['Masa', '≈ 100 millones de veces la del Sol, según Kip Thorne en «The Science of Interstellar»'],
['Ubicación', 'En otra galaxia, al otro lado del agujero de gusano que aparece cerca de Saturno'],
['Papel en la historia', 'Su gravedad domina el sistema donde orbitan los planetas candidatos de las misiones Lázaro']],
science: [['Horizonte de sucesos', 'Frontera a partir de la cual nada, ni siquiera la luz, puede escapar. Por eso la región central es completamente oscura.'],
['Disco de acreción', 'Gas que gira alrededor del agujero negro y se calienta por fricción hasta brillar intensamente. Las zonas internas son las más calientes.'],
['Lente gravitacional', 'La gravedad curva la trayectoria de la luz. Por eso vemos la parte trasera del disco por encima y por debajo del agujero negro.'],
['Dilatación temporal gravitacional', 'Cerca de una gran masa el tiempo transcurre más despacio que lejos de ella. Es un efecto real que corrigen, por ejemplo, los satélites GPS.']],
render: 'Representación en tiempo real: SOLARIS traza rayos de luz curvados por la gravedad de un agujero negro sin rotación (aproximación de Schwarzschild). Las imágenes de la película se calcularon con relatividad general para un agujero negro en rotación.',
} },
{ id: 'miller', system: 'gargantua', name: 'Planeta de Miller', short: 'Miller', type: 'planet', fiction: true, parent: 'gargantua', visR: 2.4, color: '#8fa3b3', rotH: 30,
orbit: { t: 'cine', a: 60 * 15, P: 3, phase: 0.18, i: 0.05 },
vis: { style: 'rock', feat: 'icy', c1: '#6f8596', c2: '#4b5f70', c3: '#a9b8c4', crater: 0, seed: 3.1, atm: { color: [0.62, 0.72, 0.82], h: 0.035, k: 0.9 } },
src: 'Película Interstellar (2014).',
info: {
desc: 'Mundo oceánico de aguas poco profundas, en una órbita muy próxima a Gargantúa.',
fiction: [['Tipo', 'Planeta oceánico'], ['Condiciones generales', 'Un océano poco profundo cubre casi toda la superficie; la gravedad de Gargantúa levanta olas gigantescas'],
['Relación con Gargantúa', 'Orbita muy cerca del agujero negro'], ['Dilatación temporal', 'Aproximadamente 1 hora en el planeta equivale a 7 años lejos de Gargantúa'],
['Misión', 'Primera parada de la Endurance, en busca de la Dra. Laura Miller, de las misiones Lázaro'],
['Masa, radio y temperatura', 'No especificado']],
film: 'En la película, la visita dura unas pocas horas: al regresar a la Endurance, Romilly les dice que han pasado 23 años, 4 meses y 8 días.',
dilation: { localH: 1, externalY: 7 },
} },
{ id: 'mann', system: 'gargantua', name: 'Planeta de Mann', short: 'Mann', type: 'planet', fiction: true, parent: 'gargantua', visR: 3.0, color: '#d6e0e8', rotH: 40,
orbit: { t: 'cine', a: 60 * 46, P: 40, phase: 0.62, i: -0.04 },
vis: { style: 'rock', feat: 'icy', c1: '#e3eaf0', c2: '#9eb0bf', c3: '#ffffff', crater: 0.12, seed: 8.4, atm: { color: [0.82, 0.88, 0.94], h: 0.025, k: 0.6 } },
src: 'Película Interstellar (2014).',
info: {
desc: 'Mundo helado e inhóspito, con grandes formaciones de hielo y nubes congeladas.',
fiction: [['Tipo', 'Planeta helado'], ['Condiciones generales', 'Superficie extremadamente fría, terreno irregular, niebla y nubes heladas; atmósfera no respirable'],
['Relación con Gargantúa', 'Orbita el agujero negro a mayor distancia que Miller'],
['Misión', 'Explorado por el Dr. Mann, de las misiones Lázaro, que enviaba datos alentadores sobre su habitabilidad; la tripulación descubre que no eran ciertos'],
['Masa, radio y temperatura', 'No especificado']],
} },
{ id: 'edmunds', system: 'gargantua', name: 'Planeta de Edmunds', short: 'Edmunds', type: 'planet', fiction: true, parent: 'gargantua', visR: 3.3, color: '#c9a77c', rotH: 26,
orbit: { t: 'cine', a: 60 * 80, P: 110, phase: 0.05, i: 0.03 },
vis: { style: 'rock', feat: 'marte', c1: '#a9875c', c2: '#6d5238', c3: '#d2b48a', crater: 0.04, seed: 5.6, atm: { color: [0.96, 0.84, 0.66], h: 0.03, k: 0.75 } },
src: 'Película Interstellar (2014).',
info: {
desc: 'Mundo rocoso con atmósfera que la película presenta como el candidato más favorable para la vida humana.',
fiction: [['Tipo', 'Planeta rocoso'], ['Condiciones generales', 'Terreno rocoso y abierto, con relieve y atmósfera; la película sugiere condiciones más favorables que las de Miller y Mann'],
['Relación con Gargantúa', 'Orbita el agujero negro a gran distancia'],
['Misión', 'Explorado por Wolf Edmunds, de las misiones Lázaro. Al final de la película, Amelia Brand llega allí y establece un campamento'],
['Masa, radio y temperatura', 'No especificado']],
} },
];
TYPE_LABEL.blackhole = 'Agujero negro supermasivo';
const OBLIQ = 23.4392911 * DEG;
const GAUSS_K_DEG = 0.9856076686; // movimiento medio (°/día) para a = 1 UA
const Astro = {
jdFromDate(d) { return d.getTime() / 86400000 + 2440587.5; },
dateFromJD(jd) { return new Date((jd - 2440587.5) * 86400000); },
jdFromISO(s) { const [y, m, d] = s.split('-').map(Number); return Date.UTC(y, m - 1, d, 0, 0, 0) / 86400000 + 2440587.5; },
eqToEcl(v) { const c = Math.cos(OBLIQ), s = Math.sin(OBLIQ); return [v[0], v[1] * c + v[2] * s, -v[1] * s + v[2] * c]; },
radecVec(raDeg, decDeg) { const a = raDeg * DEG, d = decDeg * DEG; return [Math.cos(d) * Math.cos(a), Math.cos(d) * Math.sin(a), Math.sin(d)]; },
radecToScene(raDeg, decDeg) { return Astro.eclToScene(Astro.eqToEcl(Astro.radecVec(raDeg, decDeg))); },
eclToScene(v) { return [v[0], v[2], -v[1]]; },
sceneToEcl(v) { return [v[0], -v[2], v[1]]; },
keplerE(M, e) {
M = wrapPi(M);
let E = e < 0.8 ? M : (M >= 0 ? Math.PI : -Math.PI);
for (let k = 0; k < 50; k++) {
const f = E - e * Math.sin(E) - M, fp = 1 - e * Math.cos(E);
const dE = f / fp; E -= dE;
if (Math.abs(dE) < 1e-12) break;
}
return E;
},
elemToEcl(a, e, iDeg, nodeDeg, periDeg, E) {
const xp = a * (Math.cos(E) - e), yp = a * Math.sqrt(1 - e * e) * Math.sin(E);
const O = nodeDeg * DEG, w = periDeg * DEG, I = iDeg * DEG;
const cO = Math.cos(O), sO = Math.sin(O), cw = Math.cos(w), sw = Math.sin(w), cI = Math.cos(I), sI = Math.sin(I);
return [
(cw * cO - sw * sO * cI) * xp + (-sw * cO - cw * sO * cI) * yp,
(cw * sO + sw * cO * cI) * xp + (-sw * sO + cw * cO * cI) * yp,
(sw * sI) * xp + (cw * sI) * yp,
];
},
elementsAt(orbit, jd) {
if (orbit.t === 'jpl') {
const T = (jd - J2000) / 36525, el = orbit.el, rt = orbit.rt;
const a = el[0] + rt[0] * T, e = el[1] + rt[1] * T, I = el[2] + rt[2] * T, L = el[3] + rt[3] * T, vp = el[4] + rt[4] * T, Om = el[5] + rt[5] * T;
return { a, e, i: I, node: Om, peri: vp - Om, M: wrap360(L - vp) };
}
if (orbit.t === 'kep') {
const n = GAUSS_K_DEG / Math.pow(orbit.a, 1.5);
return { a: orbit.a, e: orbit.e, i: orbit.i, node: orbit.node, peri: orbit.peri, M: wrap360(n * (jd - orbit.tpJD)) };
}
return null;
},
helioPos(orbit, jd) {
const k = Astro.elementsAt(orbit, jd);
const E = Astro.keplerE(k.M * DEG, k.e);
return Astro.elemToEcl(k.a, k.e, k.i, k.node, k.peri, E);
},
sampleOrbit(orbit, jd, n) {
const k = Astro.elementsAt(orbit, jd);
const pts = [], fr = [];
for (let j = 0; j <= n; j++) {
let E;
if (k.e > 0.5) {
const nu = -Math.PI + TAU * j / n;
E = 2 * Math.atan2(Math.sqrt(1 - k.e) * Math.sin(nu / 2), Math.sqrt(1 + k.e) * Math.cos(nu / 2));
} else E = -Math.PI + TAU * j / n;
const M = E - k.e * Math.sin(E);
pts.push(Astro.elemToEcl(k.a, k.e, k.i, k.node, k.peri, E));
fr.push(((M / TAU) % 1 + 1) % 1);
}
return { pts, fr, curFrac: k.M / 360 };
},
lunaGeo(jd) {
const d = jd - J2000;
const Om = wrap360(125.044 - 0.0529538083 * d), L = wrap360(218.316 + 13.176396 * d), M = wrap360(134.963 + 13.064993 * d);
const lam = L + 6.289 * Math.sin(M * DEG);
const u = (lam - Om) * DEG, O = Om * DEG, I = 5.145 * DEG;
const r = 385001 - 20905 * Math.cos(M * DEG);
return {
pos: [r * (Math.cos(O) * Math.cos(u) - Math.sin(O) * Math.sin(u) * Math.cos(I)),
r * (Math.sin(O) * Math.cos(u) + Math.cos(O) * Math.sin(u) * Math.cos(I)),
r * Math.sin(u) * Math.sin(I)],
node: Om, frac: wrap360(lam - Om) / 360,
};
},
poleFrame(raDeg, decDeg) {
const P = V.norm(Astro.radecToScene(raDeg, decDeg));
const qa = (raDeg + 90) * DEG;
const Q = V.norm(Astro.eclToScene(Astro.eqToEcl([Math.cos(qa), Math.sin(qa), 0])));
return { P, Q, PxQ: V.cross(P, Q) };
},
basisFromFrame(fr, Wdeg) {
const w = Wdeg * DEG, c = Math.cos(w), s = Math.sin(w);
const X = [fr.Q[0] * c + fr.PxQ[0] * s, fr.Q[1] * c + fr.PxQ[1] * s, fr.Q[2] * c + fr.PxQ[2] * s];
const Y = fr.P, Z = V.cross(X, Y);
return [X, Y, Z];
},
orbitFrame(Pref, Nref, iDeg, nodeDeg) {
const q = V.cross(Pref, Nref);
const O = nodeDeg * DEG, I = iDeg * DEG;
const nodeDir = V.add(V.scale(Nref, Math.cos(O)), V.scale(q, Math.sin(O)));
const perp = V.add(V.scale(Nref, -Math.sin(O)), V.scale(q, Math.cos(O)));
const ninety = V.norm(V.add(V.scale(perp, Math.cos(I)), V.scale(Pref, Math.sin(I))));
const normal = V.norm(V.cross(nodeDir, ninety));
return { u: nodeDir, v: ninety, n: normal };
},
};
const SCALES = {
visual: {
label: 'Visual', badge: 'Escala visual: tamaños ampliados y distancias comprimidas',
dist: r => 60 * Math.pow(Math.max(r, 1e-9), 0.55),
inv: s => Math.pow(Math.max(s, 0) / 60, 1 / 0.55),
body: (R, isSun) => (isSun ? 6 : Math.max(0.012, 0.55 * Math.pow(R / 6371, 0.6))),
moon: (dKm, Rp, RpS) => RpS * (1.25 + 0.75 * Math.sqrt(dKm / Rp)),
glsl: [1, 60, 0.55, 0], overview: 560, maxDist: 12000, shadowSun: 0.35,
},
real: {
label: 'Real', badge: 'Escala real: distancias y tamaños proporcionales',
dist: r => 1000 * r,
inv: s => s / 1000,
body: R => (R / AU_KM) * 1000,
moon: dKm => (dKm / AU_KM) * 1000,
glsl: [1, 1000, 1, 0], overview: 75000, maxDist: 2.5e6, shadowSun: 1,
},
edu: {
label: 'Educativa', badge: 'Escala educativa: separación casi uniforme y tamaños ampliados',
dist: r => 48 * Math.log(1 + Math.max(r, 0) / 0.15),
inv: s => 0.15 * (Math.exp(Math.max(s, 0) / 48) - 1),
body: (R, isSun) => (isSun ? 9 : Math.max(0.05, 1.1 * Math.pow(R / 6371, 0.42))),
moon: (dKm, Rp, RpS) => RpS * (1.3 + 0.55 * Math.sqrt(dKm / Rp)),
glsl: [2, 48, 0, 0.15], overview: 620, maxDist: 9000, shadowSun: 0.35,
},
};
const ScaleState = {
from: 'visual', to: 'visual', t: 1,
set(name, animate) {
if (name === this.to && this.t >= 1) return;
this.from = animate ? this.current() : name; this.to = name; this.t = animate ? 0 : 1;
if (!animate) this.from = name;
},
current() { return this.t >= 1 ? this.to : this.to; },
k() { return easeInOut(clamp(this.t, 0, 1)); },
blend(fa, fb) { const k = this.k(); if (k >= 1 || this.from === this.to) return fb; if (k <= 0) return fa; return Math.exp(lerp(Math.log(Math.max(fa, 1e-30)), Math.log(Math.max(fb, 1e-30)), k)); },
dist(r) { const A = SCALES[this.from], B = SCALES[this.to]; return this.blend(A.dist(r), B.dist(r)); },
inv(s) { return SCALES[this.to].inv(s); },
body(R, isSun) { const A = SCALES[this.from], B = SCALES[this.to]; return this.blend(A.body(R, isSun), B.body(R, isSun)); },
moon(dKm, Rp, RpSa, RpSb) { const A = SCALES[this.from], B = SCALES[this.to]; return this.blend(A.moon(dKm, Rp, RpSa), B.moon(dKm, Rp, RpSb)); },
bodyA(R, isSun) { return SCALES[this.from].body(R, isSun); },
bodyB(R, isSun) { return SCALES[this.to].body(R, isSun); },
get(key) { const A = SCALES[this.from][key], B = SCALES[this.to][key]; return this.blend(A, B); },
mapVec(v) { const r = V.len(v); if (r < 1e-12) return [0, 0, 0]; const s = this.dist(r) / r; return [v[0] * s, v[1] * s, v[2] * s]; },
};
const GLX = {
gl: null, canvas: null, W: 1, H: 1, dpr: 1, hdr: false, samples: 4,
progs: {}, scratch16: new Float32Array(16), scratch9: new Float32Array(9),
fsTri: null, targets: null, bloomLevels: 5, quality: 'alta',
init(canvas) {
const gl = canvas.getContext('webgl2', { antialias: false, alpha: false, depth: false, stencil: false, premultipliedAlpha: false, powerPreference: 'high-performance', preserveDrawingBuffer: false });
if (!gl) throw new Error('WebGL2 no disponible');
this.gl = gl; this.canvas = canvas;
this.hdr = !!gl.getExtension('EXT_color_buffer_float');
gl.getExtension('OES_texture_float_linear');
this.maxSamples = gl.getParameter(gl.MAX_SAMPLES) || 0;
this.fsTri = gl.createVertexArray(); // triángulo de pantalla completa (usa gl_VertexID)
return gl;
},
compile(type, src, name) {
const gl = this.gl, s = gl.createShader(type);
gl.shaderSource(s, src); gl.compileShader(s);
if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
const log = gl.getShaderInfoLog(s);
const lines = src.split('\n').map((l, i) => (i + 1) + ': ' + l).join('\n');
console.error('Error de shader [' + name + ']\n' + log + '\n' + lines.slice(0, 20000));
throw new Error('Shader ' + name + ': ' + log);
}
return s;
},
program(name, vs, fs) {
if (this.progs[name]) return this.progs[name];
const gl = this.gl, p = gl.createProgram();
gl.attachShader(p, this.compile(gl.VERTEX_SHADER, vs, name + ':vs'));
gl.attachShader(p, this.compile(gl.FRAGMENT_SHADER, fs, name + ':fs'));
gl.linkProgram(p);
if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error('Link ' + name + ': ' + gl.getProgramInfoLog(p));
const u = {}, n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
for (let i = 0; i < n; i++) {
const info = gl.getActiveUniform(p, i);
const base = info.name.replace(/\[0\]$/, '');
u[base] = { loc: gl.getUniformLocation(p, info.name), type: info.type, size: info.size };
}
const a = {}, na = gl.getProgramParameter(p, gl.ACTIVE_ATTRIBUTES);
for (let i = 0; i < na; i++) { const info = gl.getActiveAttrib(p, i); a[info.name] = gl.getAttribLocation(p, info.name); }
const prog = { name, p, u, a, tex: 0 };
this.progs[name] = prog;
return prog;
},
use(prog) { this.gl.useProgram(prog.p); prog.tex = 0; },
set(prog, name, v) {
const gl = this.gl, U = prog.u[name];
if (!U || v == null) return;
const L = U.loc;
switch (U.type) {
case gl.FLOAT: if (U.size > 1) gl.uniform1fv(L, v); else gl.uniform1f(L, v); break;
case gl.FLOAT_VEC2: gl.uniform2fv(L, v); break;
case gl.FLOAT_VEC3: gl.uniform3fv(L, v); break;
case gl.FLOAT_VEC4: gl.uniform4fv(L, v); break;
case gl.INT: case gl.BOOL: gl.uniform1i(L, v); break;
case gl.FLOAT_MAT4: { const s = this.scratch16; for (let i = 0; i < 16; i++) s[i] = v[i]; gl.uniformMatrix4fv(L, false, s); break; }
case gl.FLOAT_MAT3: { const s = this.scratch9; for (let i = 0; i < 9; i++) s[i] = v[i]; gl.uniformMatrix3fv(L, false, s); break; }
case gl.SAMPLER_2D: case gl.SAMPLER_CUBE: {
const unit = prog.tex++;
gl.activeTexture(gl.TEXTURE0 + unit);
gl.bindTexture(U.type === gl.SAMPLER_CUBE ? gl.TEXTURE_CUBE_MAP : gl.TEXTURE_2D, v);
gl.uniform1i(L, unit); break;
}
default: break;
}
},
setAll(prog, obj) { for (const k in obj) this.set(prog, k, obj[k]); },
mesh(attribs, index, mode, dynamic) {
const gl = this.gl, vao = gl.createVertexArray(), bufs = {};
gl.bindVertexArray(vao);
let count = 0;
for (const name in attribs) {
const at = attribs[name], b = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, b);
gl.bufferData(gl.ARRAY_BUFFER, at.data, dynamic ? gl.DYNAMIC_DRAW : gl.STATIC_DRAW);
bufs[name] = { b, size: at.size, loc: -1 };
count = at.data.length / at.size;
}
let ib = null;
if (index) { ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, index, gl.STATIC_DRAW); count = index.length; }
gl.bindVertexArray(null);
return { vao, bufs, ib, count, mode: mode == null ? gl.TRIANGLES : mode, indexType: index instanceof Uint32Array ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, bound: null };
},
freeMesh(m) { const gl = this.gl; if (!m || !m.vao) return; gl.deleteVertexArray(m.vao); for (const k in m.bufs) gl.deleteBuffer(m.bufs[k].b); if (m.ib) gl.deleteBuffer(m.ib); m.vao = null; m.count = 0; },
bindMesh(mesh, prog) {
const gl = this.gl;
gl.bindVertexArray(mesh.vao);
if (mesh.bound === prog.name) return;
for (const name in mesh.bufs) {
const B = mesh.bufs[name], loc = prog.a[name];
if (B.loc >= 0 && B.loc !== loc) gl.disableVertexAttribArray(B.loc);
if (loc == null || loc < 0) { B.loc = -1; continue; }
gl.bindBuffer(gl.ARRAY_BUFFER, B.b);
gl.enableVertexAttribArray(loc);
gl.vertexAttribPointer(loc, B.size, gl.FLOAT, false, 0, 0);
B.loc = loc;
}
if (mesh.ib) gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, mesh.ib);
mesh.bound = prog.name;
},
update(mesh, name, data, count) {
const gl = this.gl, B = mesh.bufs[name];
gl.bindBuffer(gl.ARRAY_BUFFER, B.b);
gl.bufferSubData(gl.ARRAY_BUFFER, 0, data);
if (count != null) mesh.count = count;
},
draw(mesh, prog, count) {
const gl = this.gl;
this.bindMesh(mesh, prog);
const n = count != null ? count : mesh.count;
if (mesh.ib) gl.drawElements(mesh.mode, n, mesh.indexType, 0); else gl.drawArrays(mesh.mode, 0, n);
},
state(s) {
const gl = this.gl;
if (s.depthTest) gl.enable(gl.DEPTH_TEST); else gl.disable(gl.DEPTH_TEST);
gl.depthMask(!!s.depthWrite);
if (s.cull === 'none') gl.disable(gl.CULL_FACE); else { gl.enable(gl.CULL_FACE); gl.cullFace(s.cull === 'front' ? gl.FRONT : gl.BACK); }
if (!s.blend || s.blend === 'none') gl.disable(gl.BLEND);
else {
gl.enable(gl.BLEND);
if (s.blend === 'add') gl.blendFunc(gl.ONE, gl.ONE);
else if (s.blend === 'premul') gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
else gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
}
},
texFromCanvas(cv, opts) {
const gl = this.gl, t = gl.createTexture();
gl.bindTexture(gl.TEXTURE_2D, t);
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, cv);
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, (opts && opts.clampS) ? gl.CLAMP_TO_EDGE : gl.REPEAT);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
return t;
},
mipTex: [],
avifOK() {
if (!this._avif) this._avif = new Promise(res => { const i = new Image(); i.onload = () => res(i.width === 2); i.onerror = () => res(false); i.src = 'data:image/avif;base64,AAAAIGZ0eXBhdmlmAAAAAGF2aWZtaWYxbWlhZk1BMUIAAADrbWV0YQAAAAAAAAAhaGRscgAAAAAAAAAAcGljdAAAAAAAAAAAAAAAAAAAAAAOcGl0bQAAAAAAAQAAAB5pbG9jAAAAAEQAAAEAAQAAAAEAAAETAAAAHAAAAChpaW5mAAAAAAABAAAAGmluZmUCAAAAAAEAAGF2MDFDb2xvcgAAAABqaXBycAAAAEtpcGNvAAAAFGlzcGUAAAAAAAAAAgAAAAIAAAAQcGl4aQAAAAADCAgIAAAADGF2MUOBAAwAAAAAE2NvbHJuY2x4AAEADQAGgAAAABdpcG1hAAAAAAAAAAEAAQQBAoMEAAAAJG1kYXQSAAoFGAA2BCAyERgACiiihAAAsBNU7TCk6wz8'; });
return this._avif;
},
setTexQuality(q) {
const gl = this.gl; if (!gl) return;
const base = { low: 2, medium: 1, high: 0, ultra: 0 }[q] || 0, an = gl.getExtension('EXT_texture_filter_anisotropic');
const maxA = an ? gl.getParameter(an.MAX_TEXTURE_MAX_ANISOTROPY_EXT) : 1, aniso = Math.min(maxA, { low: 1, medium: 2, high: 4, ultra: 16 }[q] || 4);
this.texQ = q;
for (const t of this.mipTex) { gl.bindTexture(gl.TEXTURE_2D, t); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_BASE_LEVEL, base); if (an) gl.texParameterf(gl.TEXTURE_2D, an.TEXTURE_MAX_ANISOTROPY_EXT, aniso); }
},
async texFromURI(uri, clampT) {
const gl = this.gl, img = new Image();
img.src = uri; await img.decode();
const t = gl.createTexture(); this.mipTex.push(t);
gl.bindTexture(gl.TEXTURE_2D, t);
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, img);
gl.generateMipmap(gl.TEXTURE_2D);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, clampT ? gl.CLAMP_TO_EDGE : gl.REPEAT);
const an = gl.getExtension('EXT_texture_filter_anisotropic');
if (an) gl.texParameterf(gl.TEXTURE_2D, an.TEXTURE_MAX_ANISOTROPY_EXT, Math.min(4, gl.getParameter(an.MAX_TEXTURE_MAX_ANISOTROPY_EXT)));
if (this.texQ) { const q = this.texQ; this.texQ = null; this.setTexQuality(q); }
return t;
},
whiteTex() {
if (this._white) return this._white;
const gl = this.gl, t = gl.createTexture();
gl.bindTexture(gl.TEXTURE_2D, t);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([255, 255, 255, 255]));
return (this._white = t);
},
drawRange(mesh, prog, first, count) {
const gl = this.gl;
this.bindMesh(mesh, prog);
gl.drawElements(mesh.mode, count, mesh.indexType, first * (mesh.indexType === gl.UNSIGNED_INT ? 4 : 2));
},
makeTex(w, h, hdr) {
const gl = this.gl, t = gl.createTexture();
gl.bindTexture(gl.TEXTURE_2D, t);
if (hdr) gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA16F, w, h, 0, gl.RGBA, gl.HALF_FLOAT, null);
else gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, w, h, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
return t;
},
fbo(tex) {
const gl = this.gl, f = gl.createFramebuffer();
gl.bindFramebuffer(gl.FRAMEBUFFER, f);
gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
return f;
},
resize(w, h, dpr) {
const gl = this.gl;
this.dpr = dpr;
const W = Math.max(2, Math.floor(w * dpr)), H = Math.max(2, Math.floor(h * dpr));
if (W === this.W && H === this.H && this.targets) return;
this.W = W; this.H = H;
this.canvas.width = W; this.canvas.height = H;
this.buildTargets();
},
buildTargets() {
const gl = this.gl, W = this.W, H = this.H, T = this.targets;
if (T) { // liberar
[T.msFbo, T.resFbo].forEach(f => f && gl.deleteFramebuffer(f));
[T.msColor, T.msDepth, T.depthRb].forEach(r => r && gl.deleteRenderbuffer(r));
gl.deleteTexture(T.sceneTex);
T.bloom.forEach(l => { gl.deleteTexture(l.tex); gl.deleteFramebuffer(l.fbo); });
}
const hdr = this.hdr, internal = hdr ? gl.RGBA16F : gl.RGBA8;
const samples = Math.min(this.samples, this.maxSamples);
const t = { samples };
t.sceneTex = this.makeTex(W, H, hdr);
t.resFbo = this.fbo(t.sceneTex);
if (samples > 0) {
t.msFbo = gl.createFramebuffer(); gl.bindFramebuffer(gl.FRAMEBUFFER, t.msFbo);
t.msColor = gl.createRenderbuffer(); gl.bindRenderbuffer(gl.RENDERBUFFER, t.msColor);
gl.renderbufferStorageMultisample(gl.RENDERBUFFER, samples, internal, W, H);
gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.RENDERBUFFER, t.msColor);
t.msDepth = gl.createRenderbuffer(); gl.bindRenderbuffer(gl.RENDERBUFFER, t.msDepth);
gl.renderbufferStorageMultisample(gl.RENDERBUFFER, samples, gl.DEPTH_COMPONENT24, W, H);
gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, t.msDepth);
if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) { t.samples = 0; }
}
if (!t.samples) {
gl.bindFramebuffer(gl.FRAMEBUFFER, t.resFbo);
t.depthRb = gl.createRenderbuffer(); gl.bindRenderbuffer(gl.RENDERBUFFER, t.depthRb);
gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT24, W, H);
gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, t.depthRb);
}
t.bloom = [];
let bw = W, bh = H;
for (let i = 0; i < this.bloomLevels; i++) {
bw = Math.max(1, bw >> 1); bh = Math.max(1, bh >> 1);
const tex = this.makeTex(bw, bh, hdr);
t.bloom.push({ tex, fbo: this.fbo(tex), w: bw, h: bh });
}
gl.bindFramebuffer(gl.FRAMEBUFFER, null);
this.targets = t;
},
beginScene() {
const gl = this.gl, t = this.targets;
gl.bindFramebuffer(gl.FRAMEBUFFER, t.samples ? t.msFbo : t.resFbo);
gl.viewport(0, 0, this.W, this.H);
gl.clearColor(0, 0, 0, 1); gl.clearDepth(1);
gl.depthMask(true);
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
},
endScene(post) {
const gl = this.gl, t = this.targets;
if (t.samples) {
gl.bindFramebuffer(gl.READ_FRAMEBUFFER, t.msFbo);
gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, t.resFbo);
gl.blitFramebuffer(0, 0, this.W, this.H, 0, 0, this.W, this.H, gl.COLOR_BUFFER_BIT, gl.NEAREST);
}
this.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'none' });
gl.bindVertexArray(this.fsTri);
const P = post.progs;
let src = t.sceneTex, sw = this.W, sh = this.H;
if (post.bloom > 0) {
for (let i = 0; i < t.bloom.length; i++) {
const L = t.bloom[i];
gl.bindFramebuffer(gl.FRAMEBUFFER, L.fbo); gl.viewport(0, 0, L.w, L.h);
const pr = i === 0 ? P.downT : P.down;
this.use(pr);
this.set(pr, 'u_src', src); this.set(pr, 'u_texel', [1 / sw, 1 / sh]); this.set(pr, 'u_threshold', post.threshold);
gl.drawArrays(gl.TRIANGLES, 0, 3);
src = L.tex; sw = L.w; sh = L.h;
}
this.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'add' });
this.use(P.up);
for (let i = t.bloom.length - 1; i > 0; i--) {
const S = t.bloom[i], D = t.bloom[i - 1];
gl.bindFramebuffer(gl.FRAMEBUFFER, D.fbo); gl.viewport(0, 0, D.w, D.h);
this.set(P.up, 'u_src', S.tex); this.set(P.up, 'u_texel', [1 / S.w, 1 / S.h]);
P.up.tex = 0;
gl.drawArrays(gl.TRIANGLES, 0, 3);
}
this.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'none' });
}
gl.bindFramebuffer(gl.FRAMEBUFFER, null);
gl.viewport(0, 0, this.W, this.H);
const C = P.comp; this.use(C);
this.set(C, 'u_scene', t.sceneTex); this.set(C, 'u_bloom', t.bloom[0].tex);
this.set(C, 'u_bloomK', post.bloom); this.set(C, 'u_exposure', post.exposure);
this.set(C, 'u_inScale', this.hdr ? 1 : 4); this.set(C, 'u_res', [this.W, this.H]); this.set(C, 'u_vignette', post.vignette);
gl.drawArrays(gl.TRIANGLES, 0, 3);
},
};
const Geo = {
sphere(nLon, nLat) {
const pos = [], idx = [];
for (let i = 0; i <= nLat; i++) {
const lat = Math.PI / 2 - Math.PI * i / nLat;
for (let j = 0; j <= nLon; j++) {
const lon = TAU * j / nLon - Math.PI;
pos.push(Math.cos(lat) * Math.cos(lon), Math.sin(lat), -Math.cos(lat) * Math.sin(lon));
}
}
const row = nLon + 1;
for (let i = 0; i < nLat; i++) for (let j = 0; j < nLon; j++) {
const a = i * row + j, b = a + 1, c = a + row, d = c + 1;
Geo.pushTri(idx, pos, a, c, b); Geo.pushTri(idx, pos, b, c, d);
}
const P = new Float32Array(pos);
return { pos: P, nrm: P, dir: P, idx: P.length / 3 > 65535 ? new Uint32Array(idx) : new Uint16Array(idx) };
},
pushTri(idx, pos, a, b, c) {
const ax = pos[a * 3], ay = pos[a * 3 + 1], az = pos[a * 3 + 2];
const e1 = [pos[b * 3] - ax, pos[b * 3 + 1] - ay, pos[b * 3 + 2] - az], e2 = [pos[c * 3] - ax, pos[c * 3 + 1] - ay, pos[c * 3 + 2] - az];
const n = V.cross(e1, e2);
const cx = ax + pos[b * 3] + pos[c * 3], cy = ay + pos[b * 3 + 1] + pos[c * 3 + 1], cz = az + pos[b * 3 + 2] + pos[c * 3 + 2];
if (n[0] * cx + n[1] * cy + n[2] * cz >= 0) idx.push(a, b, c); else idx.push(a, c, b);
},
icosphere(detail) {
const t = (1 + Math.sqrt(5)) / 2;
let verts = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]].map(V.norm);
let faces = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
for (let d = 0; d < detail; d++) {
const cache = {}, nf = [];
const mid = (a, b) => { const k = a < b ? a + '_' + b : b + '_' + a; if (cache[k] != null) return cache[k]; verts.push(V.norm(V.scale(V.add(verts[a], verts[b]), 0.5))); return (cache[k] = verts.length - 1); };
for (const [a, b, c] of faces) { const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a); nf.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]); }
faces = nf;
}
return { verts, faces };
},
irregular(shape, rough, seed, bilobe) {
const { verts, faces } = Geo.icosphere(5);
const rng = makeRng(Math.floor(seed * 1000) + 7);
const off = [rng() * 100, rng() * 100, rng() * 100];
const sx = shape ? shape[0] : 1, sy = shape ? shape[1] : 1, sz = shape ? shape[2] : 1;
const pos = new Float32Array(verts.length * 3), dir = new Float32Array(verts.length * 3);
verts.forEach((d, i) => {
let r;
if (bilobe) {
const lobe = (c, rad) => { const b = V.dot(d, c), disc = b * b - V.dot(c, c) + rad * rad; return disc >= 0 ? b + Math.sqrt(disc) : 0; };
r = Math.max(lobe([-0.42, 0, 0], 0.62), lobe([0.52, 0.04, 0], 0.5));
} else r = 1;
const n1 = Noise3.fbm(d[0] * 1.6 + off[0], d[1] * 1.6 + off[1], d[2] * 1.6 + off[2], 4);
const n2 = Noise3.fbm(d[0] * 5 + off[1], d[1] * 5 + off[2], d[2] * 5 + off[0], 3);
r *= 1 + rough * n1 + rough * 0.25 * n2;
pos[i * 3] = d[0] * r * sx; pos[i * 3 + 1] = d[1] * r * sy; pos[i * 3 + 2] = d[2] * r * sz;
dir[i * 3] = d[0]; dir[i * 3 + 1] = d[1]; dir[i * 3 + 2] = d[2];
});
let mean = 0; for (let i = 0; i < verts.length; i++) mean += Math.hypot(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
mean /= verts.length; for (let i = 0; i < pos.length; i++) pos[i] /= mean;
const idx = [];
faces.forEach(([a, b, c]) => Geo.pushTri(idx, dir, a, b, c));
const nrm = new Float32Array(pos.length);
for (let f = 0; f < idx.length; f += 3) {
const a = idx[f], b = idx[f + 1], c = idx[f + 2];
const pa = [pos[a * 3], pos[a * 3 + 1], pos[a * 3 + 2]], pb = [pos[b * 3], pos[b * 3 + 1], pos[b * 3 + 2]], pc = [pos[c * 3], pos[c * 3 + 1], pos[c * 3 + 2]];
const n = V.cross(V.sub(pb, pa), V.sub(pc, pa));
for (const k of [a, b, c]) { nrm[k * 3] += n[0]; nrm[k * 3 + 1] += n[1]; nrm[k * 3 + 2] += n[2]; }
}
for (let i = 0; i < verts.length; i++) { const l = Math.hypot(nrm[i * 3], nrm[i * 3 + 1], nrm[i * 3 + 2]) || 1; nrm[i * 3] /= l; nrm[i * 3 + 1] /= l; nrm[i * 3 + 2] /= l; }
return { pos, nrm, dir, idx: new Uint32Array(idx) };
},
ring(inner, outer, seg, rings) {
const pos = [], idx = [];
for (let r = 0; r <= rings; r++) {
const rad = inner + (outer - inner) * r / rings;
for (let s = 0; s <= seg; s++) { const a = TAU * s / seg; pos.push(Math.cos(a) * rad, 0, -Math.sin(a) * rad); }
}
const row = seg + 1;
for (let r = 0; r < rings; r++) for (let s = 0; s < seg; s++) { const a = r * row + s, b = a + 1, c = a + row, d = c + 1; idx.push(a, b, c, b, d, c); }
return { pos: new Float32Array(pos), idx: new Uint16Array(idx) };
},
quad() { return new Float32Array([-1, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1, 0]); },
};
const Noise3 = (() => {
const p = new Uint8Array(512), rng = makeRng(1337), perm = [];
for (let i = 0; i < 256; i++) perm.push(i);
for (let i = 255; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [perm[i], perm[j]] = [perm[j], perm[i]]; }
for (let i = 0; i < 512; i++) p[i] = perm[i & 255];
const fade = t => t * t * t * (t * (t * 6 - 15) + 10);
const grad = (h, x, y, z) => { const u = (h & 15) < 8 ? x : y, v = (h & 15) < 4 ? y : ((h & 15) === 12 || (h & 15) === 14) ? x : z; return ((h & 1) ? -u : u) + ((h & 2) ? -v : v); };
function noise(x, y, z) {
const X = Math.floor(x) & 255, Y = Math.floor(y) & 255, Z = Math.floor(z) & 255;
x -= Math.floor(x); y -= Math.floor(y); z -= Math.floor(z);
const u = fade(x), v = fade(y), w = fade(z);
const A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z, B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;
return lerp(lerp(lerp(grad(p[AA], x, y, z), grad(p[BA], x - 1, y, z), u), lerp(grad(p[AB], x, y - 1, z), grad(p[BB], x - 1, y - 1, z), u), v),
lerp(lerp(grad(p[AA + 1], x, y, z - 1), grad(p[BA + 1], x - 1, y, z - 1), u), lerp(grad(p[AB + 1], x, y - 1, z - 1), grad(p[BB + 1], x - 1, y - 1, z - 1), u), v), w);
}
function fbm(x, y, z, oct) { let s = 0, a = 0.5, f = 1; for (let i = 0; i < oct; i++) { s += a * noise(x * f, y * f, z * f); a *= 0.5; f *= 2.03; } return s; }
return { noise, fbm };
})();
const SH = {};
SH.HEAD = `#version 300 es
precision highp float;
precision highp int;
uniform mat4 u_view; uniform mat4 u_proj;
uniform float u_depthK; uniform float u_depthInv; uniform float u_time; uniform vec3 u_sun; uniform float u_outScale;
#define PI 3.14159265359
#define TAU 6.28318530718
#define DEG 0.01745329252
float sq(float x){ return x * x; }
`;
SH.FRAG = `
out vec4 o_col;
void writeDepth(float d){ gl_FragDepth = clamp(log2(1.0 + max(d, 0.0) * u_depthK) * u_depthInv, 0.0, 1.0); }
vec4 outc(vec3 c, float a){ return vec4(c * u_outScale, a); }
`;
SH.NOISE = `
vec3 mod289(vec3 x){ return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec4 mod289(vec4 x){ return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec4 permute(vec4 x){ return mod289(((x * 34.0) + 1.0) * x); }
vec4 taylorInvSqrt(vec4 r){ return 1.79284291400159 - 0.85373472095314 * r; }
float snoise(vec3 v){
const vec2 C = vec2(1.0/6.0, 1.0/3.0);
const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);
vec3 i = floor(v + dot(v, C.yyy));
vec3 x0 = v - i + dot(i, C.xxx);
vec3 g = step(x0.yzx, x0.xyz);
vec3 l = 1.0 - g;
vec3 i1 = min(g.xyz, l.zxy);
vec3 i2 = max(g.xyz, l.zxy);
vec3 x1 = x0 - i1 + C.xxx;
vec3 x2 = x0 - i2 + C.yyy;
vec3 x3 = x0 - D.yyy;
i = mod289(i);
vec4 p = permute(permute(permute(i.z + vec4(0.0, i1.z, i2.z, 1.0)) + i.y + vec4(0.0, i1.y, i2.y, 1.0)) + i.x + vec4(0.0, i1.x, i2.x, 1.0));
float n_ = 0.142857142857;
vec3 ns = n_ * D.wyz - D.xzx;
vec4 j = p - 49.0 * floor(p * ns.z * ns.z);
vec4 x_ = floor(j * ns.z);
vec4 y_ = floor(j - 7.0 * x_);
vec4 x = x_ * ns.x + ns.yyyy;
vec4 y = y_ * ns.x + ns.yyyy;
vec4 h = 1.0 - abs(x) - abs(y);
vec4 b0 = vec4(x.xy, y.xy);
vec4 b1 = vec4(x.zw, y.zw);
vec4 s0 = floor(b0) * 2.0 + 1.0;
vec4 s1 = floor(b1) * 2.0 + 1.0;
vec4 sh = -step(h, vec4(0.0));
vec4 a0 = b0.xzyw + s0.xzyw * sh.xxyy;
vec4 a1 = b1.xzyw + s1.xzyw * sh.zzww;
vec3 p0 = vec3(a0.xy, h.x);
vec3 p1 = vec3(a0.zw, h.y);
vec3 p2 = vec3(a1.xy, h.z);
vec3 p3 = vec3(a1.zw, h.w);
vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));
p0 *= norm.x; p1 *= norm.y; p2 *= norm.z; p3 *= norm.w;
vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), 0.0);
m = m * m;
return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}
float fbm3(vec3 p){ float s = 0.0, a = 0.5; for (int i = 0; i < 3; i++){ s += a * snoise(p); p *= 2.03; a *= 0.5; } return s; }
float fbm4(vec3 p){ float s = 0.0, a = 0.5; for (int i = 0; i < 4; i++){ s += a * snoise(p); p *= 2.03; a *= 0.5; } return s; }
float fbm5(vec3 p){ float s = 0.0, a = 0.5; for (int i = 0; i < 5; i++){ s += a * snoise(p); p *= 2.03; a *= 0.5; } return s; }
vec3 hash33(vec3 p){ p = fract(p * vec3(0.1031, 0.1030, 0.0973)); p += dot(p, p.yxz + 33.33); return fract((p.xxy + p.yxx) * p.zyx); }
vec2 worley(vec3 x){
vec3 ip = floor(x), fp = fract(x); float best = 9.0, id = 0.0;
for (int k = -1; k <= 1; k++) for (int j = -1; j <= 1; j++) for (int i = -1; i <= 1; i++){
vec3 g = vec3(float(i), float(j), float(k)); vec3 o = hash33(ip + g);
float d = length(g + o - fp); if (d < best){ best = d; id = o.z; }
}
return vec2(best, id);
}
vec3 ll(float la, float lo){ la *= DEG; lo *= DEG; return vec3(cos(la) * cos(lo), sin(la), -cos(la) * sin(lo)); }
float angDist(vec3 a, vec3 b){ return acos(clamp(dot(a, b), -1.0, 1.0)); }
float bnd(float x, float lo, float hi, float s){ return smoothstep(lo - s, lo + s, x) * (1.0 - smoothstep(hi - s, hi + s, x)); }
mat2 rot2(float a){ float c = cos(a), s = sin(a); return mat2(c, s, -s, c); }
`;
SH.VS_BODY = SH.HEAD + `
in vec3 a_pos; in vec3 a_nrm; in vec3 a_dir;
uniform mat4 u_model;
out vec3 v_obj; out vec3 v_wn; out vec3 v_wp; out float v_depth;
void main(){
vec4 wp = u_model * vec4(a_pos, 1.0);
v_wp = wp.xyz;
v_wn = normalize(mat3(u_model) * a_nrm);
v_obj = a_dir;
vec4 vp = u_view * wp;
v_depth = -vp.z;
gl_Position = u_proj * vp;
}`;
SH.LIGHT = `
in vec3 v_obj; in vec3 v_wn; in vec3 v_wp; in float v_depth;
uniform vec4 u_occ[4]; uniform int u_occN; uniform float u_sunR; uniform float u_ambient; uniform float u_hover;
uniform float u_radius; uniform vec4 u_atm; uniform float u_bump; uniform vec3 u_center;
float eclipse(vec3 wp){
float f = 1.0;
vec3 ts = u_sun - wp; float ds = length(ts); vec3 s = ts / ds;
float as = u_sunR / ds;
for (int i = 0; i < 4; i++){
if (i >= u_occN) break;
vec3 c = u_occ[i].xyz - wp; float r = u_occ[i].w;
float t = dot(c, s);
if (t <= 0.0 || t >= ds) continue;
float dp = length(c - s * t);
float ao = r / t, beta = dp / t;
float cover = 1.0 - smoothstep(abs(ao - as), ao + as, beta);
float maxD = clamp((ao * ao) / (as * as), 0.0, 1.0);
f *= 1.0 - cover * maxD;
}
return f;
}
vec3 bumpN(vec3 N, float hh){
vec3 dpx = dFdx(v_wp), dpy = dFdy(v_wp);
float dhx = dFdx(hh), dhy = dFdy(hh);
vec3 r1 = cross(dpy, N), r2 = cross(N, dpx);
float det = dot(dpx, r1);
vec3 grad = sign(det) * (dhx * r1 + dhy * r2);
return normalize(abs(det) * N - grad);
}
`;
SH.CRATER = `
float craterField(vec3 x, float seed, out float rim, out float fresh){
vec3 ip = floor(x), fp = fract(x);
float h = 0.0; rim = 0.0; fresh = 0.0;
for (int k = -1; k <= 1; k++) for (int j = -1; j <= 1; j++) for (int i = -1; i <= 1; i++){
vec3 g = vec3(float(i), float(j), float(k));
vec3 o = hash33(ip + g + seed);
if (o.z > 0.62) continue;
float rad = 0.16 + 0.34 * o.x * o.x;
vec3 c = g + 0.15 + 0.7 * o - fp;
float d = length(c) / rad;
if (d > 2.6) continue;
float bowl = d < 1.0 ? (d * d - 1.0) * 0.85 : 0.0;
float r1 = exp(-sq((d - 1.0) * 3.6));
float ej = d > 1.0 ? exp(-(d - 1.0) * 2.4) * 0.12 : 0.0;
float peak = rad > 0.33 ? exp(-d * d * 45.0) * 0.4 : 0.0;
h += (bowl + r1 * 0.32 + ej + peak) * rad;
rim = max(rim, r1);
fresh = max(fresh, step(0.86, o.y) * exp(-max(d - 0.9, 0.0) * 1.6));
}
return h;
}
`;
SH.FS_ROCK = (feat) => SH.HEAD + '#define F_' + feat.toUpperCase() + '\n' + SH.FRAG + SH.NOISE + SH.LIGHT + SH.CRATER + `
uniform vec3 u_c1; uniform vec3 u_c2; uniform vec3 u_c3; uniform float u_crater; uniform float u_seed; uniform float u_detail;
void main(){
vec3 p = normalize(v_obj);
float lat = asin(clamp(p.y, -1.0, 1.0));
float lon = atan(-p.z, p.x);
float n = fbm4(p * 2.5 + u_seed);
float n2 = fbm4(p * 9.0 + u_seed * 1.37);
vec3 alb = mix(u_c2, u_c1, smoothstep(-0.45, 0.45, n));
alb = mix(alb, u_c3, smoothstep(0.25, 0.75, n2) * 0.35);
float h = n * 0.22 + n2 * 0.07;
float rim1, fr1, rim2, fr2, rim3, fr3;
float crA = craterField(p * 3.0, u_seed, rim1, fr1);
float crB = craterField(p * 8.0, u_seed + 17.0, rim2, fr2);
float crC = 0.0; rim3 = 0.0; fr3 = 0.0;
if (u_detail > 0.5) crC = craterField(p * 22.0, u_seed + 31.0, rim3, fr3);
float cr = crA + crB * 0.45 + crC * 0.18;
float fresh = max(fr1, max(fr2 * 0.7, fr3 * 0.45)) * u_crater;
float craterK = u_crater;
float spec = 0.0; float lunar = 0.75; float emit = 0.0;
#ifdef F_LUNA
vec3 MAR[16] = vec3[16](vec3(32.8,-15.6,0.17), vec3(28.0,17.5,0.11), vec3(8.5,31.4,0.13), vec3(17.0,59.1,0.085), vec3(-7.8,51.3,0.1), vec3(-15.2,35.5,0.06), vec3(-21.3,-16.6,0.1), vec3(22.0,-55.0,0.2), vec3(2.0,-52.0,0.17), vec3(38.0,-45.0,0.1), vec3(56.0,-15.0,0.07), vec3(56.0,15.0,0.065), vec3(-24.4,-38.6,0.07), vec3(7.5,-30.9,0.08), vec3(13.3,3.6,0.05), vec3(-10.0,-23.0,0.06));
float mare = 0.0;
for (int i = 0; i < 16; i++){ float d = angDist(p, ll(MAR[i].x, MAR[i].y)) + 0.06 * n2 + 0.04 * n; mare = max(mare, 1.0 - smoothstep(MAR[i].z * 0.7, MAR[i].z * 1.08, d)); }
alb = mix(alb, u_c2 * (0.9 + 0.2 * n2), mare * 0.9);
craterK *= 1.0 - mare * 0.65;
vec3 RC[2] = vec3[2](vec3(-43.3, -11.2, 0.55), vec3(9.6, -20.1, 0.3));
for (int i = 0; i < 2; i++){
vec3 c = ll(RC[i].x, RC[i].y); float d = angDist(p, c);
vec3 t1 = normalize(cross(c, vec3(0.0, 1.0, 0.0))); vec3 t2 = cross(c, t1);
vec3 v = p - c * dot(p, c); float ang = atan(dot(v, t2), dot(v, t1));
float ray = smoothstep(0.35, 0.85, snoise(vec3(cos(ang) * 9.0, sin(ang) * 9.0, d * 2.0 + float(i) * 7.0)));
float ry = ray * exp(-d / RC[i].z) * smoothstep(0.015, 0.05, d) + (1.0 - smoothstep(0.012, 0.025, d));
alb = mix(alb, u_c3 * 1.08, clamp(ry, 0.0, 1.0) * 0.65);
}
#endif
#ifdef F_MERCURIO
alb = mix(alb, u_c3 * 1.05, fresh * 0.6);
alb *= 0.92 + 0.12 * smoothstep(-0.3, 0.4, fbm4(p * 4.0 + 9.0));
#endif
#ifdef F_MARTE
float dark = smoothstep(0.05, 0.35, fbm4(p * 2.2 + 5.0) + 0.25 * fbm4(p * 7.0 + 2.0));
dark = max(dark, 1.0 - smoothstep(0.12, 0.34, angDist(p, ll(8.4, 69.5)) + 0.12 * n2));
dark = max(dark, (1.0 - smoothstep(0.2, 0.45, abs(lat + 0.35) + 0.15 * n)) * 0.6);
float hellas = 1.0 - smoothstep(0.12, 0.3, angDist(p, ll(-42.4, 70.5)) + 0.05 * n);
alb = mix(u_c1 * (0.92 + 0.16 * n2), u_c2, clamp(dark, 0.0, 1.0) * 0.72);
alb = mix(alb, u_c3, hellas * 0.55);
float vmLat = -8.0 * DEG + 0.035 * sin(lon * 3.0);
float inLon = smoothstep(-118.0 * DEG, -102.0 * DEG, lon) * (1.0 - smoothstep(-45.0 * DEG, -30.0 * DEG, lon));
float vm = exp(-sq((lat - vmLat) / 0.022)) * inLon;
h -= vm * 0.7; alb = mix(alb, u_c2 * 0.7, vm * 0.5);
float om = angDist(p, ll(18.65, -133.8)); h += exp(-om * om / 0.0035) * 0.9;
vec3 TH[3] = vec3[3](vec3(-8.4, -120.5, 0.0), vec3(1.6, -112.6, 0.0), vec3(11.8, -104.5, 0.0));
for (int i = 0; i < 3; i++){ float dv = angDist(p, ll(TH[i].x, TH[i].y)); h += exp(-dv * dv / 0.0016) * 0.6; }
float cap = smoothstep(1.27, 1.35, lat + 0.05 * n2) + smoothstep(1.33, 1.40, -lat + 0.05 * n2);
alb = mix(alb, vec3(0.92, 0.9, 0.88), clamp(cap, 0.0, 1.0));
craterK *= 1.0 - clamp(cap, 0.0, 1.0);
#endif
#ifdef F_IO
alb = mix(u_c1, vec3(0.78, 0.55, 0.2), smoothstep(0.0, 0.6, n) * 0.6);
alb = mix(alb, u_c3, smoothstep(0.3, 0.7, n2) * 0.55);
vec2 w = worley(p * 7.0 + u_seed);
float sel = fract(w.y * 13.7);
alb = mix(alb, vec3(0.05, 0.04, 0.03), 1.0 - smoothstep(0.03, 0.075, w.x));
alb = mix(alb, vec3(0.62, 0.2, 0.08), (1.0 - smoothstep(0.16, 0.34, w.x)) * smoothstep(0.06, 0.12, w.x) * step(0.55, sel));
alb *= mix(1.0, 0.62, smoothstep(0.85, 1.3, abs(lat)));
lunar = 0.4;
#endif
#ifdef F_EUROPA
float l1 = pow(max(1.0 - abs(snoise(p * 3.0 + 2.0)), 0.0), 18.0);
float l2 = pow(max(1.0 - abs(snoise(p * 7.0 + 9.0)), 0.0), 26.0);
float l3 = pow(max(1.0 - abs(snoise(p * 14.0 + 4.0)), 0.0), 34.0);
float lin = clamp(l1 + l2 * 0.8 + l3 * 0.6, 0.0, 1.0);
alb = mix(u_c1, u_c3, smoothstep(-0.3, 0.4, n) * 0.4);
alb = mix(alb, mix(u_c2, u_c1, 0.45), smoothstep(0.3, 0.6, n2) * 0.45);
alb = mix(alb, u_c2, lin * 0.75);
h += lin * 0.05; lunar = 0.4;
#endif
#ifdef F_GANIMEDES
float terr = smoothstep(-0.1, 0.2, fbm4(p * 2.0 + 3.0));
float grooves = 0.5 + 0.5 * sin((p.x * 0.6 + p.y * 0.3 + p.z * 0.7) * 95.0 + n2 * 22.0);
alb = mix(u_c2, u_c3 * (0.86 + 0.14 * grooves), terr);
alb = mix(alb, u_c3 * 1.1, fresh * 0.7);
craterK *= 1.0 - terr * 0.4;
#endif
#ifdef F_CALISTO
alb = mix(u_c2, u_c1, smoothstep(-0.4, 0.4, n));
alb = mix(alb, u_c3, clamp(fresh + rim2 * 0.35 + rim3 * 0.25, 0.0, 1.0) * 0.75);
float dv = angDist(p, ll(15.0, -56.0));
alb = mix(alb, u_c3 * 0.8, (1.0 - smoothstep(0.04, 0.16, dv)) * 0.7 + (0.5 + 0.5 * cos(dv * 70.0)) * 0.18 * (1.0 - smoothstep(0.25, 0.85, dv)));
#endif
#ifdef F_ENCELADO
float south = 1.0 - smoothstep(-1.28, -1.05, lat);
vec2 q = p.xz; float s = dot(q, normalize(vec2(1.0, 0.6))) + 0.03 * snoise(p * 6.0);
float fs = abs(fract(s * 9.0) - 0.5);
float line = (1.0 - smoothstep(0.02, 0.06, fs)) * (1.0 - smoothstep(0.18, 0.25, abs(s)));
alb = mix(alb, vec3(0.4, 0.55, 0.65), line * south * 0.75);
craterK *= 1.0 - smoothstep(-0.9, -0.3, -lat) * 0.0 - (1.0 - smoothstep(-0.9, -0.2, lat)) * 0.8;
lunar = 0.3;
#endif
#ifdef F_JAPETO
float dk = 1.0 - smoothstep(1.0, 1.35, angDist(p, ll(0.0, -90.0)) + 0.25 * n);
dk *= 1.0 - smoothstep(0.8, 1.2, abs(lat) + 0.1 * n2);
alb = mix(u_c1, u_c2, dk);
h += exp(-sq(lat / 0.035)) * 0.5 * (0.6 + 0.4 * n2);
#endif
#ifdef F_TRITON
float cap = 1.0 - smoothstep(-0.35, -0.05, lat + 0.15 * n);
vec2 wc = worley(p * 14.0);
alb = mix(mix(u_c2, u_c1, 0.5 + 0.3 * wc.x), u_c3, cap);
h += (1.0 - cap) * wc.x * 0.25;
float st = smoothstep(0.6, 0.85, snoise(vec3(p.x * 22.0, p.y * 3.0, p.z * 22.0)));
alb = mix(alb, u_c2 * 0.6, st * cap * 0.5);
lunar = 0.4;
#endif
#ifdef F_PLUTON
float lonE = lon < 0.0 ? lon + TAU : lon;
float spn = 1.0 - smoothstep(0.22, 0.33, angDist(p, ll(25.0, 177.0)) + 0.05 * n2);
float east = 1.0 - smoothstep(0.22, 0.4, angDist(p, ll(-5.0, 215.0)) + 0.08 * n);
float heart = max(spn, east * 0.8);
float cth = (1.0 - smoothstep(0.28, 0.55, abs(lat + 0.12) + 0.12 * n)) * smoothstep(10.0 * DEG, 35.0 * DEG, lonE) * (1.0 - smoothstep(145.0 * DEG, 168.0 * DEG, lonE));
alb = mix(u_c1 * (0.9 + 0.2 * n2), u_c2, cth * 0.85);
alb = mix(alb, u_c3, heart);
vec2 wp2 = worley(p * 30.0);
alb *= mix(1.0, 0.9 + 0.1 * smoothstep(0.0, 0.1, wp2.x), spn);
craterK *= 1.0 - heart;
lunar = 0.35;
#endif
#ifdef F_CARONTE
float capc = smoothstep(0.95, 1.25, lat + 0.12 * n);
alb = mix(alb, u_c2, capc * 0.8);
#endif
#ifdef F_CERES
float oc = angDist(p, ll(19.8, -120.7));
alb = mix(alb, vec3(0.92), 1.0 - smoothstep(0.008, 0.03, oc + 0.01 * n2));
alb = mix(alb, vec3(0.8), (1.0 - smoothstep(0.004, 0.015, angDist(p, ll(20.6, -119.0)))) * 0.8);
#endif
#ifdef F_VESTA
float rs = angDist(p, ll(-75.0, -59.0)) / 0.9;
h += (rs < 1.0 ? (rs * rs - 1.0) * 0.8 : 0.0) + exp(-rs * rs * 30.0) * 0.6 + exp(-sq((rs - 1.0) * 5.0)) * 0.2;
alb = mix(alb, u_c3, smoothstep(0.0, 0.6, fbm4(p * 3.0 + 4.0)) * 0.3);
#endif
#ifdef F_METAL
spec = 0.25;
#endif
#ifdef F_RUBBLE
vec2 wr = worley(p * 26.0);
h += (1.0 - smoothstep(0.0, 0.28, wr.x)) * 0.35;
alb *= 0.85 + 0.3 * wr.y;
#endif
#ifdef F_MIMAS
float hd = angDist(p, ll(1.66, -111.8)) / 0.33;
h += (hd < 1.0 ? (hd * hd - 1.0) * 0.75 : 0.0) + exp(-sq((hd - 1.0) * 4.0)) * 0.3 + exp(-hd * hd * 60.0) * 0.55;
#endif
#ifdef F_MIRANDA
vec3 CO[3] = vec3[3](vec3(-60.0, 30.0, 0.0), vec3(-40.0, 180.0, 0.0), vec3(-20.0, -60.0, 0.0));
for (int i = 0; i < 3; i++){
float dc = angDist(p, ll(CO[i].x, CO[i].y));
float m = 1.0 - smoothstep(0.45, 0.6, dc + 0.08 * n);
float stripes = 0.5 + 0.5 * cos(dc * 70.0 + n2 * 6.0);
alb = mix(alb, mix(u_c2, u_c3, stripes), m * 0.85);
craterK *= 1.0 - m * 0.7;
}
#endif
#ifdef F_TITAN
alb = mix(u_c1, u_c2, smoothstep(0.1, 0.5, fbm4(p * 2.0 + 1.0)) * 0.22);
alb = mix(alb, u_c3, smoothstep(0.9, 1.4, abs(lat)) * 0.25);
craterK = 0.0; lunar = 0.0;
#endif
#ifdef F_COMET
lunar = 0.5;
#endif
h += cr * craterK;
if (craterK > 0.0) alb = mix(alb, u_c3, fresh * 0.4);
vec3 Ng = normalize(v_wn);
float fpw = length(fwidth(p));
float bk = u_bump * u_radius * (1.0 - smoothstep(0.006, 0.06, fpw));
vec3 N = bumpN(Ng, h * bk);
vec3 L = normalize(u_sun - v_wp); vec3 Vv = normalize(-v_wp);
float ndl = dot(N, L), ndv = max(dot(N, Vv), 0.0);
float lam = max(ndl, 0.0);
float ls = ndl > 0.0 ? 2.0 * ndl / (ndl + ndv + 0.05) : 0.0;
float diff = mix(lam, min(ls, 1.4) * 0.72, lunar);
diff *= smoothstep(-0.08, 0.12, dot(Ng, L));
float sh = eclipse(v_wp);
vec3 col = alb * diff * sh * 1.35 + alb * u_ambient;
if (spec > 0.0){ vec3 hv = normalize(L + Vv); col += spec * pow(max(dot(N, hv), 0.0), 40.0) * sh * step(0.0, ndl); }
float fres = pow(max(1.0 - max(dot(Ng, Vv), 0.0), 0.0), 3.0);
float litA = smoothstep(-0.3, 0.35, dot(Ng, L));
col += u_atm.rgb * fres * litA * u_atm.a;
col += vec3(0.55, 0.7, 1.0) * fres * u_hover * 0.35;
writeDepth(v_depth);
o_col = outc(col, 1.0);
}`;
SH.FS_EARTH = SH.HEAD + SH.FRAG + SH.NOISE + SH.LIGHT + `
uniform sampler2D u_land;
void main(){
vec3 p = normalize(v_obj);
float lat = asin(clamp(p.y, -1.0, 1.0)), lon = atan(-p.z, p.x);
float nA = fbm4(p * 6.0), nB = fbm4(p * 6.0 + 7.3);
vec2 uv = vec2(lon / TAU + 0.5, lat / PI + 0.5) + vec2(nA, nB) * 0.0035;
vec3 tx = texture(u_land, uv).rgb;
float nC = fbm5(p * 38.0);
float land = smoothstep(0.47, 0.53, tx.r + nC * 0.12);
float inland = tx.g, mount = tx.b;
float ad = abs(lat) / DEG;
float desert = smoothstep(12.0, 20.0, ad + nB * 6.0) * (1.0 - smoothstep(30.0, 38.0, ad + nA * 6.0)) * smoothstep(0.38, 0.68, inland + nA * 0.35);
float rain = 1.0 - smoothstep(6.0, 20.0, ad + nA * 10.0);
float tundra = smoothstep(55.0, 66.0, ad + nB * 8.0);
float lonD = lon / DEG;
float green = step(0.0, lat) * smoothstep(58.0, 62.0, ad) * smoothstep(-58.0, -52.0, lonD) * (1.0 - smoothstep(-25.0, -18.0, lonD));
float ice = smoothstep(70.0, 75.0, ad + nB * 5.0) + green + (1.0 - step(0.0, lat)) * smoothstep(60.0, 63.0, ad);
vec3 cForest = vec3(0.035, 0.07, 0.025), cGrass = vec3(0.09, 0.11, 0.045), cDesert = vec3(0.42, 0.29, 0.15), cTundra = vec3(0.13, 0.12, 0.09), cRock = vec3(0.17, 0.13, 0.1), cIce = vec3(0.8, 0.84, 0.9);
vec3 lc = mix(cGrass, cForest, rain);
lc = mix(lc, cGrass * 1.15, smoothstep(0.2, 0.7, fbm4(p * 15.0)) * 0.5);
lc = mix(lc, cDesert * (0.9 + 0.25 * nC), clamp(desert, 0.0, 1.0));
lc = mix(lc, cTundra, tundra);
lc = mix(lc, cRock, smoothstep(0.2, 0.8, mount) * 0.55);
float snow = clamp(ice + smoothstep(0.78, 1.0, mount + nC * 0.25) * smoothstep(25.0, 45.0, ad) * 0.7, 0.0, 1.0);
lc = mix(lc, cIce, snow);
vec3 deep = vec3(0.003, 0.014, 0.05), shallow = vec3(0.008, 0.05, 0.1);
vec3 ocn = mix(deep, shallow, smoothstep(0.2, 0.5, inland));
float seaIce = clamp(step(0.0, lat) * smoothstep(77.0, 82.0, ad + nB * 6.0) + (1.0 - step(0.0, lat)) * smoothstep(64.0, 70.0, ad + nB * 5.0), 0.0, 1.0);
ocn = mix(ocn, vec3(0.7, 0.76, 0.82), seaIce);
vec3 alb = mix(ocn, lc, land);
float ocean = (1.0 - land) * (1.0 - seaIce);
float h = land * (mount * 0.6 + nC * 0.12);
vec3 Ng = normalize(v_wn);
float fpw = length(fwidth(p));
vec3 N = bumpN(Ng, h * u_bump * u_radius * (1.0 - smoothstep(0.004, 0.04, fpw)));
vec3 L = normalize(u_sun - v_wp), Vv = normalize(-v_wp);
float ndl = dot(N, L), ndv = max(dot(Ng, Vv), 0.0);
float day = smoothstep(-0.12, 0.18, dot(Ng, L));
float sh = eclipse(v_wp);
vec3 col = alb * max(ndl, 0.0) * smoothstep(-0.08, 0.12, dot(Ng, L)) * sh * 1.5 + alb * u_ambient;
vec3 hv = normalize(L + Vv);
col += vec3(1.0, 0.9, 0.75) * pow(max(dot(Ng, hv), 0.0), 260.0) * 0.9 * ocean * sh * step(0.0, ndl);
col += vec3(0.45, 0.55, 0.65) * pow(max(dot(Ng, hv), 0.0), 30.0) * 0.035 * ocean * sh * step(0.0, ndl);
float cityN = smoothstep(0.62, 0.95, fbm5(p * 55.0 + 3.0) * 0.5 + 0.5 + (1.0 - inland) * 0.22);
float city = land * (1.0 - snow) * (1.0 - clamp(desert, 0.0, 1.0) * 0.85) * cityN * (1.0 - smoothstep(55.0, 65.0, ad));
col += vec3(1.0, 0.66, 0.32) * city * (1.0 - day) * 0.5;
float fres = pow(max(1.0 - ndv, 0.0), 2.5);
col = mix(col, u_atm.rgb * 0.5 * day * sh, clamp(fres * 0.7 + 0.08 * day, 0.0, 1.0) * u_atm.a);
col += vec3(0.55, 0.7, 1.0) * fres * u_hover * 0.35;
writeDepth(v_depth);
o_col = outc(col, 1.0);
}`;
SH.FS_CLOUDS = SH.HEAD + SH.FRAG + SH.NOISE + SH.LIGHT + `
uniform float u_cloudT;
void main(){
vec3 p = normalize(v_obj);
float lat = asin(clamp(p.y, -1.0, 1.0));
vec3 q = p; q.xz = rot2(u_cloudT) * q.xz;
vec3 w = q + 0.35 * vec3(fbm3(q * 2.0), fbm3(q * 2.0 + 3.3), fbm3(q * 2.0 + 6.1));
float c = fbm5(w * 3.3);
float band = 0.17 * exp(-sq(lat / 0.12)) + 0.16 * exp(-sq((abs(lat) - 0.95) / 0.16)) - 0.16 * exp(-sq((abs(lat) - 0.42) / 0.12));
c = smoothstep(0.06, 0.55, c + band);
c *= 0.75 + 0.35 * fbm3(w * 16.0);
c = clamp(c, 0.0, 1.0);
vec3 Ng = normalize(v_wn); vec3 L = normalize(u_sun - v_wp);
float diff = smoothstep(-0.12, 0.4, dot(Ng, L));
float sh = eclipse(v_wp);
vec3 col = vec3(1.0, 0.99, 0.97) * diff * sh * 1.25 + vec3(u_ambient);
float a = c * 0.92;
writeDepth(v_depth);
o_col = outc(col * a, a);
}`;
SH.FS_VENUS = SH.HEAD + SH.FRAG + SH.NOISE + SH.LIGHT + `
void main(){
vec3 p = normalize(v_obj);
float lat = asin(clamp(p.y, -1.0, 1.0));
vec3 q = p; q.xz = rot2(u_time * 0.015) * q.xz;
float lon = atan(-q.z, q.x);
vec3 w = q * vec3(1.0, 2.6, 1.0) + 0.5 * vec3(fbm3(q * 3.0), fbm3(q * 3.0 + 2.0), 0.0);
float n = fbm5(w * 2.0);
float chev = sin(lon * 2.0 + abs(lat) * 5.0 + n * 3.0);
vec3 alb = mix(vec3(0.5, 0.36, 0.17), vec3(0.8, 0.68, 0.45), smoothstep(-0.6, 0.6, n + chev * 0.22));
alb = mix(alb, vec3(0.85, 0.79, 0.65), smoothstep(0.9, 1.35, abs(lat)) * 0.5);
vec3 Ng = normalize(v_wn); vec3 L = normalize(u_sun - v_wp); vec3 Vv = normalize(-v_wp);
float ndl = dot(Ng, L);
float diff = max(ndl * 0.85 + 0.15, 0.0) * smoothstep(-0.2, 0.15, ndl);
float sh = eclipse(v_wp);
vec3 col = alb * diff * sh * 1.35 + alb * u_ambient;
float fres = pow(max(1.0 - max(dot(Ng, Vv), 0.0), 0.0), 2.5);
col += u_atm.rgb * fres * smoothstep(-0.3, 0.3, ndl) * u_atm.a * 0.6;
col += vec3(0.55, 0.7, 1.0) * fres * u_hover * 0.35;
writeDepth(v_depth);
o_col = outc(col, 1.0);
}`;
SH.RING = `
vec4 ringSaturn(float r, float fw){
float s = max(0.0025, fw * 1.5);
float n1 = snoise(vec3(r * 170.0, 0.5, 0.5)); float n2 = snoise(vec3(r * 900.0, 1.5, 0.5));
float st = 0.82 + 0.12 * n1 + 0.06 * n2;
float D = bnd(r, 1.11, 1.236, s) * 0.03;
float C = bnd(r, 1.239, 1.527, s) * (0.1 + 0.07 * smoothstep(1.4, 1.52, r));
float B = bnd(r, 1.527, 1.951, s) * (0.62 + 0.33 * smoothstep(1.53, 1.75, r));
float CD = bnd(r, 1.951, 2.027, s) * 0.06;
float A = bnd(r, 2.027, 2.269, s) * 0.55;
float enck = bnd(r, 2.2135, 2.2195, max(0.0008, fw));
float F = bnd(r, 2.322, 2.33, max(0.0015, fw)) * 0.45;
float a = (D + C + B + CD + A * (1.0 - enck)) * st + F;
vec3 col = mix(vec3(0.42, 0.36, 0.3), vec3(0.85, 0.76, 0.6), smoothstep(1.45, 1.6, r));
col = mix(col, vec3(0.7, 0.66, 0.6), smoothstep(1.98, 2.05, r));
return vec4(col * (0.9 + 0.15 * n1), clamp(a, 0.0, 1.0));
}
vec4 ringUranus(float r, float fw){
float s = max(0.0018, fw * 1.5); float a = 0.0;
a += bnd(r, 1.635, 1.640, s); a += bnd(r, 1.650, 1.655, s); a += bnd(r, 1.664, 1.669, s);
a += bnd(r, 1.748, 1.753, s); a += bnd(r, 1.785, 1.790, s); a += bnd(r, 1.844, 1.848, s);
a += bnd(r, 1.861, 1.866, s); a += bnd(r, 1.888, 1.893, s); a += bnd(r, 1.994, 2.006, s) * 1.6;
return vec4(vec3(0.55, 0.56, 0.6), clamp(a * 0.55, 0.0, 1.0));
}
vec4 ringNeptune(float r, float fw){
float s = max(0.003, fw * 1.5);
float a = bnd(r, 1.68, 1.72, 0.01) * 0.25 + bnd(r, 2.145, 2.152, s) * 0.8 + bnd(r, 2.16, 2.30, 0.02) * 0.08 + bnd(r, 2.535, 2.545, s) * 1.0;
return vec4(vec3(0.55, 0.5, 0.48), clamp(a, 0.0, 1.0));
}
vec4 ringFaint(float r, float rin, float rout){
float a = bnd(r, rin, rout, (rout - rin) * 0.15) * (0.6 + 0.4 * snoise(vec3(r * 60.0, 0.0, 0.0)));
return vec4(vec3(0.7, 0.6, 0.5), clamp(a, 0.0, 1.0));
}
`;
SH.FS_GAS = (feat) => SH.HEAD + '#define F_' + feat.toUpperCase() + '\n' + SH.FRAG + SH.NOISE + SH.LIGHT + SH.RING + `
uniform vec3 u_ringN; uniform float u_ringOn;
float ringShadow(vec3 wp, vec3 L){
float dn = dot(L, u_ringN);
if (abs(dn) < 1e-5) return 1.0;
float t = dot(u_center - wp, u_ringN) / dn;
if (t <= 0.0) return 1.0;
float r = length(wp + L * t - u_center) / u_radius;
return 1.0 - ringSaturn(r, 0.0).a * 0.85;
}
void main(){
vec3 p = normalize(v_obj);
float lat = asin(clamp(p.y, -1.0, 1.0));
float t = u_time * 0.01;
float jet = sin(lat * 18.0) * 0.5 + sin(lat * 7.0) * 0.3;
vec3 q = p; q.xz = rot2(jet * t) * q.xz;
vec3 w = q + 0.055 * vec3(snoise(q * 5.0 + 1.0), snoise(q * 5.0 + 3.0), snoise(q * 5.0 + 5.0)) + 0.018 * vec3(snoise(q * 16.0 + 2.0), snoise(q * 16.0 + 4.0), snoise(q * 16.0 + 6.0));
float Lt = asin(clamp(w.y, -1.0, 1.0)) / DEG;
float lonq = atan(-w.z, w.x) / DEG;
float streak = snoise(vec3(w.x * 6.0, w.y * 60.0, w.z * 6.0));
vec3 c;
#ifdef F_JUPITER
vec3 zone = vec3(0.62, 0.53, 0.40), belt = vec3(0.36, 0.22, 0.12), belt2 = vec3(0.46, 0.34, 0.22), polar = vec3(0.34, 0.31, 0.28);
c = zone;
c = mix(c, belt, bnd(Lt, 7.0, 18.0, 1.6));
c = mix(c, belt * 1.05, bnd(Lt, -20.0, -7.0, 1.6));
c = mix(c, belt2, bnd(Lt, 24.0, 31.0, 2.0));
c = mix(c, belt2, bnd(Lt, -35.0, -28.0, 2.0));
c = mix(c, belt2 * 0.95, bnd(Lt, 36.0, 42.0, 2.0) * 0.7);
c = mix(c, belt2 * 0.95, bnd(Lt, -45.0, -39.0, 2.0) * 0.6);
c = mix(c, polar, smoothstep(45.0, 64.0, abs(Lt)));
c = mix(c, vec3(0.78, 0.55, 0.32), bnd(Lt, -6.0, 6.0, 3.0) * 0.35);
c *= 0.88 + 0.14 * streak;
c = mix(c, c * vec3(1.05, 0.96, 0.9), smoothstep(-0.3, 0.6, snoise(w * 3.0 + 4.0)) * 0.6);
float fest = smoothstep(0.55, 0.8, snoise(vec3(w.x * 9.0, w.y * 20.0, w.z * 9.0))) * bnd(Lt, 3.0, 8.0, 1.5);
c = mix(c, vec3(0.2, 0.2, 0.24), fest * 0.6);
float ov = smoothstep(0.62, 0.78, snoise(w * 11.0)) * bnd(Lt, -42.0, -32.0, 2.0);
c = mix(c, vec3(0.85, 0.82, 0.76), ov * 0.8);
float dlon = mod(lonq - 60.0 + 180.0, 360.0) - 180.0;
float dla = Lt + 22.0;
float e = length(vec2(dlon * cos(22.0 * DEG) / 12.0, dla / 6.5));
float sw = (1.0 - e) * 4.0;
vec2 lc = rot2(sw) * vec2(dlon, dla);
float grsT = snoise(vec3(lc * 0.35, 3.0));
float grs = 1.0 - smoothstep(0.75, 1.0, e);
c = mix(c, vec3(0.55, 0.2, 0.09) * (0.85 + 0.3 * grsT), grs);
c = mix(c, zone * 1.08, smoothstep(0.92, 1.0, e) * (1.0 - smoothstep(1.0, 1.3, e)) * 0.55);
#endif
#ifdef F_SATURNO
vec3 zone = vec3(0.68, 0.56, 0.34), belt = vec3(0.55, 0.42, 0.24), polar = vec3(0.36, 0.40, 0.40);
c = mix(zone, belt, (0.5 + 0.5 * sin(Lt * 0.32 + 0.6)) * 0.55 + bnd(Lt, 16.0, 28.0, 3.0) * 0.45 + bnd(Lt, -30.0, -18.0, 3.0) * 0.45 + bnd(Lt, 36.0, 44.0, 2.0) * 0.3);
c = mix(c, vec3(0.75, 0.65, 0.42), bnd(Lt, -10.0, 10.0, 4.0) * 0.35);
c = mix(c, polar, smoothstep(55.0, 72.0, Lt) * 0.7 + smoothstep(60.0, 78.0, -Lt) * 0.5);
c *= 0.95 + 0.06 * streak;
float colat = 90.0 - asin(clamp(p.y, -1.0, 1.0)) / DEG;
float th = atan(p.z, p.x);
float rh = cos(PI / 6.0) / cos(mod(th, PI / 3.0) - PI / 6.0);
float hex = 1.0 - smoothstep(12.0 * rh - 0.7, 12.0 * rh + 0.7, colat);
c = mix(c, vec3(0.30, 0.38, 0.40), hex * 0.75);
c = mix(c, vec3(0.12, 0.14, 0.15), 1.0 - smoothstep(0.8, 2.2, colat));
#endif
#ifdef F_URANO
c = vec3(0.36, 0.66, 0.70);
c *= 0.97 + 0.03 * sin(Lt * 0.35) + 0.015 * streak;
c = mix(c, vec3(0.52, 0.76, 0.78), smoothstep(40.0, 65.0, Lt) * 0.55);
#endif
#ifdef F_NEPTUNO
c = vec3(0.07, 0.17, 0.62);
c = mix(c, vec3(0.045, 0.11, 0.45), bnd(Lt, -55.0, -22.0, 5.0) * 0.6 + bnd(Lt, 50.0, 75.0, 6.0) * 0.4);
c *= 0.94 + 0.08 * streak;
float cl = smoothstep(0.5, 0.75, snoise(vec3(w.x * 3.0, w.y * 26.0, w.z * 3.0))) * (bnd(Lt, -48.0, -30.0, 4.0) + bnd(Lt, 22.0, 34.0, 4.0) * 0.7);
c = mix(c, vec3(0.8, 0.85, 0.95), cl * 0.8);
float dlon = mod(lonq - 100.0 + 180.0, 360.0) - 180.0;
float e = length(vec2(dlon * cos(20.0 * DEG) / 13.0, (Lt + 20.0) / 6.0));
c = mix(c, vec3(0.02, 0.05, 0.25), (1.0 - smoothstep(0.7, 1.0, e)) * 0.8);
float e2 = length(vec2((dlon + 6.0) * cos(20.0 * DEG) / 7.0, (Lt + 27.0) / 2.5));
c = mix(c, vec3(0.85, 0.88, 0.95), (1.0 - smoothstep(0.5, 1.0, e2)) * 0.7);
#endif
vec3 Ng = normalize(v_wn); vec3 L = normalize(u_sun - v_wp); vec3 Vv = normalize(-v_wp);
float ndl = dot(Ng, L), ndv = max(dot(Ng, Vv), 0.0);
float diff = max(ndl, 0.0) * smoothstep(-0.05, 0.1, ndl);
float sh = eclipse(v_wp);
if (u_ringOn > 0.5) sh *= ringShadow(v_wp, L);
vec3 col = c * diff * sh * 1.45 * (0.65 + 0.35 * pow(max(ndv, 0.0), 0.35)) + c * u_ambient;
float fres = pow(max(1.0 - ndv, 0.0), 3.0);
col += u_atm.rgb * fres * smoothstep(-0.25, 0.3, ndl) * u_atm.a;
col += vec3(0.55, 0.7, 1.0) * fres * u_hover * 0.35;
writeDepth(v_depth);
o_col = outc(col, 1.0);
}`;
SH.SPH = `
uniform sampler2D u_map;
vec4 sph(sampler2D t, vec3 p, float du){
float lat = asin(clamp(p.y, -1.0, 1.0)), lon = atan(-p.z, p.x);
vec2 uv = vec2(lon / TAU + 0.5 + du, 0.5 - lat / PI);
float uB = fract(uv.x + 0.5) - 0.5;
vec2 dxA = dFdx(uv), dyA = dFdy(uv); float dxB = dFdx(uB), dyB = dFdy(uB);
vec2 dx = vec2(abs(dxA.x) < abs(dxB) ? dxA.x : dxB, dxA.y);
vec2 dy = vec2(abs(dyA.x) < abs(dyB) ? dyA.x : dyB, dyA.y);
return textureGrad(t, uv, dx, dy);
}
`;
SH.FS_PTEX = (feat) => SH.HEAD + '#define F_' + feat.toUpperCase() + '\n' + SH.FRAG + SH.NOISE + SH.LIGHT + SH.SPH + `
uniform sampler2D u_nmap; uniform float u_nScale; uniform mat3 u_rot; uniform sampler2D u_land; uniform float u_lunar;
uniform sampler2D u_ringTex; uniform vec2 u_rr2; uniform vec3 u_ringN; uniform float u_ringOn;
float ringA(float r){ float s = (r - u_rr2.x) / (u_rr2.y - u_rr2.x); return (s < 0.0 || s > 1.0) ? 0.0 : texture(u_ringTex, vec2(s, 0.5)).a; }
void main(){
vec3 p = normalize(v_obj);
float lat = asin(clamp(p.y, -1.0, 1.0)), lon = atan(-p.z, p.x);
float du = 0.0;
#ifdef F_GAS
du = (sin(lat * 18.0) * 0.5 + sin(lat * 7.0) * 0.3) * u_time * 0.00035;
#endif
vec3 sr = sph(u_map, p, du).rgb;
vec3 alb = pow(sr, vec3(2.2));
vec3 Ng = normalize(v_wn), N = Ng;
vec3 L = normalize(u_sun - v_wp), Vv = normalize(-v_wp);
float fpw = length(fwidth(p)), lunar = u_lunar;
#ifdef F_ROCK
float hgt = dot(sr, vec3(0.3, 0.59, 0.11));
N = bumpN(Ng, hgt * u_bump * u_radius * (1.0 - smoothstep(0.003, 0.03, fpw)));
#endif
#ifdef F_LUNA
vec3 nm = sph(u_nmap, p, 0.0).xyz * 2.0 - 1.0;
vec3 T = vec3(-sin(lon), 0.0, -cos(lon)), Bt = vec3(-sin(lat) * cos(lon), cos(lat), sin(lat) * sin(lon));
vec3 no = normalize(T * nm.x * u_nScale + Bt * nm.y * u_nScale + p * max(nm.z, 0.2));
N = normalize(mix(Ng, normalize(u_rot * no), 1.0 - smoothstep(0.004, 0.04, fpw)));
#endif
#ifdef F_VENUS
alb = mix(alb, vec3(0.62, 0.48, 0.27), 0.3);
#endif
#ifdef F_TITAN
alb = mix(alb, vec3(0.58, 0.34, 0.12), 0.42);
#endif
float ndl = dot(N, L), ndv = max(dot(N, Vv), 0.0);
float sh = eclipse(v_wp);
#ifdef F_RINGTEX
if (u_ringOn > 0.5){ float dn = dot(L, u_ringN); if (abs(dn) > 1e-5){ float t = dot(u_center - v_wp, u_ringN) / dn; if (t > 0.0) sh *= 1.0 - ringA(length(v_wp + L * t - u_center) / u_radius) * 0.85; } }
#endif
float lam = max(ndl, 0.0);
float ls = ndl > 0.0 ? 2.0 * ndl / (ndl + ndv + 0.05) : 0.0;
float diff = mix(lam, min(ls, 1.4) * 0.72, lunar) * smoothstep(-0.08, 0.12, dot(Ng, L));
#ifdef F_GAS
diff = max(dot(Ng, L), 0.0) * smoothstep(-0.05, 0.1, dot(Ng, L)) * (0.68 + 0.32 * pow(max(dot(Ng, Vv), 0.0), 0.35));
#endif
vec3 col = alb * diff * sh * 1.4 + alb * u_ambient;
#ifdef F_EARTH
float ocean = smoothstep(0.15, 0.35, sr.b - sr.r) * (1.0 - smoothstep(0.62, 0.8, min(sr.r, min(sr.g, sr.b))));
vec3 hv = normalize(L + Vv);
col += vec3(1.0, 0.9, 0.75) * pow(max(dot(Ng, hv), 0.0), 260.0) * 0.9 * ocean * sh * step(0.0, ndl);
col += vec3(0.45, 0.55, 0.65) * pow(max(dot(Ng, hv), 0.0), 30.0) * 0.035 * ocean * sh * step(0.0, ndl);
vec3 lt = texture(u_land, vec2(lon / TAU + 0.5, lat / PI + 0.5)).rgb;
float cloud = smoothstep(0.55, 0.85, min(sr.r, min(sr.g, sr.b)));
float cityN = smoothstep(0.62, 0.95, fbm5(p * 55.0 + 3.0) * 0.5 + 0.5 + (1.0 - lt.g) * 0.22);
float day = smoothstep(-0.12, 0.18, dot(Ng, L));
col += vec3(1.0, 0.66, 0.32) * lt.r * (1.0 - cloud) * cityN * (1.0 - smoothstep(55.0, 65.0, abs(lat) / DEG)) * (1.0 - day) * 0.5;
#endif
float fres = pow(max(1.0 - max(dot(Ng, Vv), 0.0), 0.0), 3.0);
col += u_atm.rgb * fres * smoothstep(-0.3, 0.35, dot(Ng, L)) * u_atm.a;
col += vec3(0.55, 0.7, 1.0) * fres * u_hover * 0.35;
writeDepth(v_depth);
o_col = outc(col, 1.0);
}`;
SH.FS_SUNTEX = SH.HEAD + SH.FRAG + SH.NOISE + SH.SPH + `
in vec3 v_obj; in vec3 v_wn; in vec3 v_wp; in float v_depth;
uniform float u_sunI;
void main(){
vec3 p = normalize(v_obj);
vec3 N = normalize(v_wn), Vv = normalize(-v_wp);
float mu = clamp(dot(N, Vv), 0.0, 1.0);
vec3 base = pow(sph(u_map, p, u_time * 0.0006).rgb, vec3(2.2));
float g1 = snoise(p * 70.0 + vec3(0.0, u_time * 0.03, 0.0)), g2 = snoise(p * 150.0 - vec3(u_time * 0.05, 0.0, 0.0));
float limb = 0.35 + 0.65 * pow(mu, 0.5);
vec3 col = base * (0.86 + 0.18 * g1 + 0.1 * g2) * limb * mix(vec3(1.0, 0.62, 0.35), vec3(1.0), smoothstep(0.0, 0.9, mu));
writeDepth(v_depth);
o_col = outc(col * u_sunI * 1.6, 1.0);
}`;
SH.VS_RING = SH.HEAD + `
in vec3 a_pos; uniform mat4 u_model;
out vec3 v_wp; out float v_r; out float v_depth;
void main(){ vec4 wp = u_model * vec4(a_pos, 1.0); v_wp = wp.xyz; v_r = length(a_pos.xz); vec4 vp = u_view * wp; v_depth = -vp.z; gl_Position = u_proj * vp; }`;
SH.FS_RING = (type) => SH.HEAD + '#define R_' + type.toUpperCase() + '\n' + SH.FRAG + SH.NOISE + SH.RING + `
in vec3 v_wp; in float v_r; in float v_depth;
uniform vec3 u_center; uniform float u_radius; uniform vec3 u_ringN; uniform float u_ambient; uniform float u_alpha; uniform vec2 u_rr; uniform sampler2D u_ringTex; uniform vec2 u_rr2;
void main(){
float fw = fwidth(v_r);
#ifdef R_SATURN
vec4 pr = ringSaturn(v_r, fw);
#endif
#ifdef R_URANUS
vec4 pr = ringUranus(v_r, fw);
#endif
#ifdef R_FAINT
vec4 pr = ringFaint(v_r, u_rr.x, u_rr.y);
#endif
#ifdef R_NEPTUNE
vec4 pr = ringNeptune(v_r, fw);
#endif
#ifdef R_SATURNTEX
float sR = (v_r - u_rr2.x) / (u_rr2.y - u_rr2.x);
vec4 tR = (sR < 0.0 || sR > 1.0) ? vec4(0.0) : texture(u_ringTex, vec2(sR, 0.5));
vec4 pr = vec4(pow(tR.rgb, vec3(2.2)) * (0.92 + 0.08 * snoise(vec3(v_r * 600.0, 1.5, 0.5))), tR.a);
#endif
float a = pr.a * u_alpha;
if (a < 0.002) discard;
vec3 L = normalize(u_sun - v_wp); vec3 Vv = normalize(-v_wp);
vec3 oc = v_wp - u_center; float b = dot(oc, L);
float sh = 1.0;
if (b < 0.0){ float dp = sqrt(max(dot(oc, oc) - b * b, 0.0)); sh = smoothstep(u_radius * 0.985, u_radius * 1.02, dp); }
float sl = dot(u_ringN, L), vl = dot(u_ringN, Vv);
bool same = sl * vl > 0.0;
float lit = same ? 1.0 : (0.28 * (1.0 - pr.a) + 0.05);
float elev = smoothstep(0.0, 0.25, abs(sl)) * 0.6 + 0.4;
float phase = same ? 1.0 : 1.0 + 1.5 * pow(max(dot(-L, Vv), 0.0), 6.0);
vec3 col = pr.rgb * (lit * elev * sh * phase * 1.25 + u_ambient);
writeDepth(v_depth);
o_col = outc(col * a, a);
}`;
SH.FS_ATM = SH.HEAD + SH.FRAG + `
in vec3 v_obj; in vec3 v_wn; in vec3 v_wp; in float v_depth;
uniform vec3 u_center; uniform float u_R; uniform float u_H; uniform vec4 u_atm; uniform float u_k;
void main(){
vec3 rd = normalize(v_wp);
float tc = dot(u_center, rd);
vec3 clos = rd * tc - u_center; float hh = length(clos);
float x = (hh - u_R) / u_H;
float top = u_R + u_H;
float chord = 2.0 * sqrt(max(top * top - hh * hh, 0.0)) / top;
float dens = exp(-max(x, 0.0) * 4.0 * u_k) * (1.0 - smoothstep(0.8, 1.0, x));
if (hh < u_R && tc > 0.0) dens *= 0.0;
vec3 n = normalize(clos + 1e-6);
vec3 Ls = normalize(u_sun - u_center);
float cl = dot(n, Ls);
float lit = smoothstep(-0.38, 0.22, cl);
float mu = dot(rd, Ls);
float phase = 0.75 * (1.0 + mu * mu) + 0.8 * pow(max(mu, 0.0), 10.0);
vec3 col = u_atm.rgb * dens * chord * lit * phase * u_atm.a * 1.6;
float term = smoothstep(-0.38, -0.05, cl) * (1.0 - smoothstep(-0.05, 0.25, cl));
col = mix(col, col * vec3(1.6, 0.65, 0.35), term * 0.7);
writeDepth(v_depth);
o_col = outc(col, 0.0);
}`;
SH.FS_SUN = SH.HEAD + SH.FRAG + SH.NOISE + `
in vec3 v_obj; in vec3 v_wn; in vec3 v_wp; in float v_depth;
uniform float u_sunI;
void main(){
vec3 p = normalize(v_obj);
float t = u_time;
vec3 N = normalize(v_wn); vec3 Vv = normalize(-v_wp);
float mu = clamp(dot(N, Vv), 0.0, 1.0);
float g1 = snoise(p * 55.0 + vec3(0.0, t * 0.03, 0.0));
float g2 = snoise(p * 130.0 - vec3(t * 0.05, 0.0, 0.0));
float gran = 0.5 + 0.42 * g1 + 0.22 * g2;
float sg = fbm3(p * 9.0 + t * 0.01);
float lat = asin(clamp(p.y, -1.0, 1.0));
float belt = smoothstep(0.08, 0.2, abs(lat)) * (1.0 - smoothstep(0.45, 0.65, abs(lat)));
float sp = snoise(p * 5.0 + vec3(13.1, 2.0, t * 0.002)) + 0.15 * snoise(p * 18.0);
float spot = smoothstep(0.7, 0.76, sp) * belt;
float umbra = smoothstep(0.76, 0.82, sp) * belt;
float fac = smoothstep(0.45, 0.6, sp) * (1.0 - smoothstep(0.6, 0.7, sp)) * belt * (1.0 - mu) * 2.0;
float limb = 0.35 + 0.65 * pow(max(mu, 0.0), 0.5);
vec3 hot = vec3(1.0, 0.74, 0.40), warm = vec3(0.95, 0.36, 0.08);
vec3 col = mix(warm, hot, smoothstep(0.0, 0.9, mu)) * (0.82 + 0.3 * gran + 0.1 * sg) * limb;
col *= 1.0 - spot * 0.55 - umbra * 0.35;
col += fac * vec3(1.0, 0.75, 0.45) * 0.35;
writeDepth(v_depth);
o_col = outc(col * u_sunI, 1.0);
}`;
SH.VS_BILL = SH.HEAD + `
in vec3 a_pos; uniform vec3 u_center; uniform float u_size;
out vec2 v_uv; out float v_depth;
void main(){
vec3 right = vec3(u_view[0][0], u_view[1][0], u_view[2][0]);
vec3 up = vec3(u_view[0][1], u_view[1][1], u_view[2][1]);
vec3 wp = u_center + (right * a_pos.x + up * a_pos.y) * u_size;
vec4 vp = u_view * vec4(wp, 1.0); v_depth = -vp.z; v_uv = a_pos.xy; gl_Position = u_proj * vp;
}`;
SH.FS_CORONA = SH.HEAD + SH.FRAG + SH.NOISE + `
in vec2 v_uv; in float v_depth;
uniform float u_scaleR; uniform float u_int;
void main(){
float q = length(v_uv);
if (q > 1.0) discard;
float r = q * u_scaleR;
float ang = atan(v_uv.y, v_uv.x);
float glow = 0.6 * exp(-(r - 1.0) * 2.4) + 0.28 * exp(-(r - 1.0) * 0.55) + 0.05 / (r * r);
float st = snoise(vec3(cos(ang) * 2.6, sin(ang) * 2.6, r * 0.35 - u_time * 0.02));
float st2 = snoise(vec3(cos(ang) * 9.0, sin(ang) * 9.0, r * 0.8 + u_time * 0.01));
glow *= 1.0 + 0.9 * sq(max(st, 0.0)) * exp(-(r - 1.0) * 0.4) + 0.25 * st2 * exp(-(r - 1.0) * 1.5);
glow *= smoothstep(0.96, 1.03, r) * (1.0 - smoothstep(0.55, 1.0, q));
vec3 col = vec3(1.0, 0.7, 0.4) * glow * u_int;
writeDepth(v_depth);
o_col = outc(col, 0.0);
}`;
SH.FS_GLOW = SH.HEAD + SH.FRAG + `
in vec2 v_uv; in float v_depth;
uniform vec3 u_color; uniform float u_int; uniform float u_sharp; uniform float u_near;
void main(){
float r = length(v_uv);
if (r > 1.0) discard;
float g = (exp(-r * r * u_sharp) + 0.25 * exp(-r * 3.5) * (1.0 - r)) * (u_near > 0.0 ? smoothstep(u_near * 0.6, u_near * 3.0, v_depth) : 1.0);
writeDepth(v_depth);
o_col = outc(u_color * g * u_int, 0.0);
}`;
SH.VS_TAIL = SH.HEAD + `
in vec3 a_pos; in vec2 a_uv; out vec2 v_uv; out float v_depth;
void main(){ vec4 vp = u_view * vec4(a_pos, 1.0); v_depth = -vp.z; v_uv = a_uv; gl_Position = u_proj * vp; }`;
SH.FS_TAIL = SH.HEAD + SH.FRAG + SH.NOISE + `
in vec2 v_uv; in float v_depth;
uniform vec3 u_color; uniform float u_int; uniform float u_seed; uniform float u_near;
void main(){
float s = clamp(v_uv.x, 0.0, 1.0), t = clamp(v_uv.y, -1.0, 1.0);
float along = pow(max(1.0 - s, 0.0), 1.6) * smoothstep(0.0, 0.03, s);
float across = exp(-t * t * 3.5);
float streak = 0.65 + 0.35 * snoise(vec3(t * 5.0, s * 3.0 - u_time * 0.05, u_seed));
float a = along * across * streak * u_int * smoothstep(u_near * 0.5, u_near * 3.0, v_depth);
writeDepth(v_depth);
o_col = outc(u_color * a, 0.0);
}`;
SH.VS_ORBIT = SH.HEAD + `
in vec3 a_pos; in float a_frac; uniform mat4 u_model;
out float v_frac; out float v_depth;
void main(){ vec4 wp = u_model * vec4(a_pos, 1.0); vec4 vp = u_view * wp; v_depth = -vp.z; v_frac = a_frac; gl_Position = u_proj * vp; }`;
SH.FS_ORBIT = SH.HEAD + SH.FRAG + `
in float v_frac; in float v_depth;
uniform vec3 u_color; uniform float u_alpha; uniform float u_cur; uniform float u_fade;
void main(){
float d = fract(u_cur - v_frac);
float tail = 0.16 + 0.84 * pow(max(1.0 - d, 0.0), 3.0);
float a = u_alpha * mix(1.0, tail, u_fade);
writeDepth(v_depth);
o_col = outc(u_color, a);
}`;
SH.VS_HZ = SH.HEAD + `
in vec2 a_pos; uniform mat4 u_model; uniform float u_rIn; uniform float u_rOut;
out float v_t; out float v_depth;
void main(){
float r = mix(u_rIn, u_rOut, a_pos.x);
vec4 wp = u_model * vec4(cos(a_pos.y) * r, 0.0, -sin(a_pos.y) * r, 1.0);
vec4 vp = u_view * wp; v_depth = -vp.z; v_t = a_pos.x; gl_Position = u_proj * vp;
}`;
SH.FS_HZ = SH.HEAD + SH.FRAG + `
in float v_t; in float v_depth;
uniform vec3 u_color; uniform vec3 u_edge; uniform float u_alpha;
void main(){
float w = max(fwidth(v_t), 1e-4);
float e = max(1.0 - smoothstep(w * 0.5, w * 2.5, v_t), smoothstep(1.0 - w * 2.5, 1.0 - w * 0.5, v_t));   // bordes nítidos de 1–2 px
float g = 0.6 + 0.4 * sin(3.14159 * v_t);                                                                // gradiente suave hacia el centro
writeDepth(v_depth);
o_col = outc(mix(u_color, u_edge, e), u_alpha * (0.2 * g + 0.65 * e));
}`;
SH.FS_BH = SH.HEAD + SH.FRAG + SH.NOISE + `
in vec2 v_uv;
uniform vec3 u_camR; uniform vec3 u_camU; uniform vec3 u_camF; uniform vec2 u_tan; uniform vec2 u_off; uniform vec3 u_cam;
uniform int u_steps; uniform float u_rin; uniform float u_rout; uniform float u_reveal; uniform float u_heat; uniform float u_spin;
float h3(vec3 p){ p = fract(p * 0.3183099 + 0.1); p *= 17.0; return fract(p.x * p.y * p.z * (p.x + p.y + p.z)); }
vec3 lensStars(vec3 d){
vec3 q = d * 220.0, c = floor(q); float h = h3(c), s = smoothstep(0.988, 1.0, h);
s *= smoothstep(0.42, 0.0, length(fract(q) - 0.5));
return s * mix(vec3(1.0, 0.86, 0.72), vec3(0.78, 0.86, 1.0), h3(c + 7.0)) * 1.6;
}
vec4 disk(vec3 p, float rr){
float edge = smoothstep(u_rin, u_rin * 1.12, rr) * (1.0 - smoothstep(u_rout * 0.62, u_rout, rr));
float a = u_spin * u_time / pow(rr, 1.5);
vec2 q = mat2(cos(a), -sin(a), sin(a), cos(a)) * p.xz;
float n1 = fbm4(vec3(q * 0.75, rr * 0.35)) * 0.5 + 0.5;
float n2 = fbm4(vec3(q * 2.6, rr * 1.7 + 3.0)) * 0.5 + 0.5;
float bands = 0.75 + 0.25 * sin(rr * 9.0 + n1 * 6.0);
float T = pow(u_rin / rr, 0.75) * pow(max(1.0 - sqrt(u_rin / rr), 0.0), 0.25) * 2.4;
float I = T * (0.45 + 0.85 * n1 * n2) * bands * edge;
vec3 col = mix(vec3(1.0, 0.42, 0.10), vec3(1.0, 0.82, 0.55), clamp(T * 1.35, 0.0, 1.0));
col = mix(col, vec3(1.0, 0.95, 0.86), clamp(T * T * 1.2 - 0.5, 0.0, 1.0));   // tono cálido: blanco solo en la zona más interna
return vec4(col * I * 1.55 * u_heat, clamp(I * 1.45, 0.0, 1.0));
}
void main(){
vec2 ndc = v_uv * 2.0 - 1.0;
vec3 dir = normalize(u_camF + (ndc.x + u_off.x) * u_tan.x * u_camR + (ndc.y + u_off.y) * u_tan.y * u_camU);   // incluye el desplazamiento óptico
vec3 x = u_cam;
float RI = u_rout * 1.7;                                   // región de influencia apreciable
float tca = -dot(x, dir), b = length(cross(x, dir));
if (b > RI || (tca < 0.0 && length(x) > RI)) discard;       // descarte barato: el rayo nunca se acerca
if (length(x) > RI) x += dir * max(tca - sqrt(max(RI * RI - b * b, 0.0)), 0.0);   // avanzar hasta la región de influencia
vec3 v = dir, hv = cross(x, v); float h2 = dot(hv, hv);
vec3 col = vec3(0.0); float alpha = 0.0; bool cap = false;
float k = 150.0 / float(u_steps);
for (int i = 0; i < 320; i++) {
if (i >= u_steps) break;
float r = length(x);
float dt = clamp(0.06 * r, 0.025, 1.4) * k;
vec3 vn = v - 1.5 * h2 * x / pow(r, 5.0) * dt;           // ecuación de la órbita de un fotón (Binet)
vec3 xn = x + vn * dt;
if (x.y * xn.y < 0.0) {                                  // cruce del plano del disco
vec3 p = mix(x, xn, x.y / (x.y - xn.y)); float rr = length(p.xz);
if (rr > u_rin && rr < u_rout) { vec4 d = disk(p, rr); col += (1.0 - alpha) * d.rgb; alpha += (1.0 - alpha) * d.a; if (alpha > 0.985) break; }
}
x = xn; v = vn;
if (dot(x, x) < 1.0) { cap = true; break; }              // horizonte de sucesos
if (dot(x, x) > RI * RI * 1.6 && dot(x, v) > 0.0) break; // escapa
}
float defl = acos(clamp(dot(normalize(v), dir), -1.0, 1.0));
vec3 bg = cap ? vec3(0.0) : lensStars(normalize(v)) * smoothstep(0.02, 0.2, defl);
float aBg = cap ? 1.0 : smoothstep(0.035, 0.45, defl) * 0.92;   // cerca de la sombra, el cielo lenteado sustituye al real
float aOut = max(alpha, aBg);
o_col = outc((col + (1.0 - alpha) * bg) * u_reveal, aOut * u_reveal);
}`;
SH.VS_STARS = SH.HEAD + `
in vec3 a_pos; in vec4 a_col; uniform float u_dpr; uniform float u_bright;
out vec3 v_col;
void main(){
vec3 vd = mat3(u_view) * a_pos;
vec4 cp = u_proj * vec4(vd, 1.0);
cp.z = cp.w * 0.99999;
gl_Position = cp;
float b = pow(10.0, -0.4 * (a_col.w - 2.0) * 0.8);
gl_PointSize = clamp(1.4 + 2.0 * sqrt(b), 1.3, 8.0) * u_dpr;
v_col = a_col.rgb * min(b, 4.0) * u_bright;
}`;
SH.FS_STARS = SH.HEAD + SH.FRAG + `
in vec3 v_col;
void main(){
vec2 d = gl_PointCoord * 2.0 - 1.0; float r2 = dot(d, d);
float g = exp(-r2 * 4.5);
if (g < 0.01) discard;
o_col = outc(v_col * g, 1.0);
}`;
SH.VS_LINESKY = SH.HEAD + `
in vec3 a_pos;
void main(){ vec3 vd = mat3(u_view) * a_pos; vec4 cp = u_proj * vec4(vd, 1.0); cp.z = cp.w * 0.99999; gl_Position = cp; }`;
SH.FS_LINESKY = SH.HEAD + SH.FRAG + `
uniform vec3 u_color; uniform float u_alpha;
void main(){ o_col = outc(u_color, u_alpha); }`;
SH.FS_SKYBAKE = `#version 300 es
precision highp float; precision highp int;
#define PI 3.14159265359
#define TAU 6.28318530718
#define DEG 0.01745329252
float sq(float x){ return x * x; }
out vec4 o_col;
uniform int u_face; uniform float u_size;
uniform vec3 u_gnp; uniform vec3 u_gc; uniform vec3 u_gy;
uniform vec4 u_pDir[12]; uniform vec4 u_pCol[12]; uniform vec4 u_pAx[12]; uniform int u_pN;
` + SH.NOISE + `
vec3 skyColor(vec3 d){
float b = asin(clamp(dot(d, u_gnp), -1.0, 1.0));
float l = atan(dot(d, u_gy), dot(d, u_gc));
float bd = b / DEG, ld = l / DEG;
float thick = 6.5 + 4.5 * (1.0 - smoothstep(0.0, 90.0, abs(ld)));
float disk = exp(-sq(bd / thick));
float halo = exp(-sq(bd / (thick * 2.6))) * 0.25;
float lumL = 0.3 + 0.7 * sq(cos(l * 0.5));
float bulge = exp(-sq(ld / 13.0) - sq(bd / 9.0));
float n1 = fbm5(d * 5.0 + 3.0); float n2 = fbm5(d * 14.0 + 9.0); float n3 = fbm4(d * 42.0 + 1.0);
float clouds = clamp(0.55 + 0.65 * n1 + 0.4 * n2 + 0.2 * n3, 0.0, 1.7);
float mw = (disk * clouds + halo) * lumL + bulge * 1.3 * (0.7 + 0.45 * n2);
float rift = exp(-sq((bd - 1.5 - 2.5 * n1) / 2.4)) * smoothstep(-40.0, -8.0, ld) * (1.0 - smoothstep(55.0, 85.0, ld));
float lanes = smoothstep(0.05, 0.45, fbm4(d * 9.0 + 1.7)) * disk;
float dust = clamp(rift * 0.85 + lanes * 0.55, 0.0, 0.95);
mw *= 1.0 - dust;
vec3 cd = mix(vec3(0.6, 0.7, 1.0), vec3(1.0, 0.84, 0.62), clamp(bulge * 1.6 + lumL * 0.35 - 0.12, 0.0, 1.0));
vec3 col = cd * mw * 0.06;
col += vec3(0.32, 0.2, 0.13) * dust * disk * 0.012;
col += vec3(0.006, 0.008, 0.014) * (0.6 + 0.4 * fbm3(d * 2.0));
for (int i = 0; i < 12; i++){
if (i >= u_pN) break;
vec3 c = u_pDir[i].xyz; float rad = u_pDir[i].w;
float cd2 = dot(d, c);
if (cd2 < cos(min(rad * 3.0, 1.5))) continue;
vec3 ax = u_pAx[i].xyz; float el = u_pAx[i].w;
vec3 off = d - c * cd2; float x = dot(off, ax); float y = length(off - ax * x);
float r = length(vec2(x / el, y)) / sin(rad);
float g = exp(-r * r * 2.2);
float tex = 0.55 + 0.7 * fbm4(d * 70.0 + float(i) * 3.0);
float inten = u_pCol[i].w;
if (inten < 0.0) col *= 1.0 - g * 0.8;
else col += u_pCol[i].rgb * g * tex * inten * 0.07;
}
vec3 gp = d * 95.0; vec3 gi = floor(gp); vec3 hh = hash33(gi + 11.0);
if (hh.x > 0.94){ vec3 cp = gi + hh; float dd = length(gp - cp); col += vec3(0.9, 0.85, 0.8) * exp(-dd * dd * 55.0) * 0.016 * hh.y; }
return col;
}
void main(){
vec2 f = gl_FragCoord.xy / u_size * 2.0 - 1.0; float u = f.x, v = f.y; vec3 d;
if (u_face == 0) d = vec3(1.0, -v, -u);
else if (u_face == 1) d = vec3(-1.0, -v, u);
else if (u_face == 2) d = vec3(u, 1.0, v);
else if (u_face == 3) d = vec3(u, -1.0, -v);
else if (u_face == 4) d = vec3(u, -v, 1.0);
else d = vec3(-u, -v, -1.0);
o_col = vec4(skyColor(normalize(d)), 1.0);
}`;
SH.VS_FULL = `#version 300 es
out vec2 v_uv;
void main(){ vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2)); v_uv = p; gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0); }`;
SH.FS_SKY = SH.HEAD + SH.FRAG + `
in vec2 v_uv;
uniform samplerCube u_sky; uniform vec3 u_camR; uniform vec3 u_camU; uniform vec3 u_camF; uniform vec2 u_tan; uniform vec2 u_off; uniform float u_skyB; uniform vec2 u_shift;
void main(){
vec2 ndc = v_uv * 2.0 - 1.0 + u_shift;
vec3 d = normalize(u_camF + u_camR * (ndc.x + u_off.x) * u_tan.x + u_camU * (ndc.y + u_off.y) * u_tan.y);
o_col = outc(texture(u_sky, d).rgb * u_skyB, 1.0);
}`;
SH.VS_BELT = SH.HEAD + `
in vec4 a_e1; in vec4 a_e2;
uniform float u_days; uniform vec4 u_mapA; uniform vec4 u_mapB; uniform float u_mapT; uniform vec3 u_off; uniform float u_dpr; uniform float u_ptK;
out vec3 v_col; out float v_a; out float v_depth;
float mapR(vec4 m, float r){ return m.x > 1.5 ? m.y * log(1.0 + r / m.w) : m.y * pow(r, m.z); }
void main(){
float a = a_e1.x, e = a_e1.y, inc = a_e1.z, node = a_e1.w, peri = a_e2.x;
float n = 0.01720209895 / (a * sqrt(a));
float M = mod(a_e2.y + n * u_days, TAU);
float E = M + e * sin(M);
E = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
E = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
E = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
float xp = a * (cos(E) - e), yp = a * sqrt(1.0 - e * e) * sin(E);
float cO = cos(node), sO = sin(node), cw = cos(peri), sw = sin(peri), cI = cos(inc), sI = sin(inc);
vec3 ecl = vec3((cw * cO - sw * sO * cI) * xp + (-sw * cO - cw * sO * cI) * yp, (cw * sO + sw * cO * cI) * xp + (-sw * sO + cw * cO * cI) * yp, sw * sI * xp + cw * sI * yp);
float r = length(ecl);
float s = exp(mix(log(mapR(u_mapA, r)), log(mapR(u_mapB, r)), u_mapT));
vec3 sc = vec3(ecl.x, ecl.z, -ecl.y) * (s / r);
vec4 vp = u_view * vec4(sc + u_off, 1.0);
gl_Position = u_proj * vp;
v_depth = -vp.z;
float px = u_ptK * a_e2.z / max(v_depth, 1e-6);
gl_PointSize = clamp(px, 1.0, 4.5) * u_dpr;
v_a = clamp(0.35 + px * 0.6, 0.35, 1.0);
float tint = a_e2.w;
v_col = tint < 0.5 ? mix(vec3(0.55, 0.5, 0.45), vec3(0.62, 0.52, 0.4), tint * 2.0) : mix(vec3(0.5, 0.55, 0.62), vec3(0.65, 0.6, 0.7), tint * 2.0 - 1.0);
}`;
SH.FS_BELT = SH.HEAD + SH.FRAG + `
in vec3 v_col; in float v_a; in float v_depth; uniform float u_alpha;
void main(){
vec2 d = gl_PointCoord * 2.0 - 1.0; float r2 = dot(d, d);
if (r2 > 1.0) discard;
float a = (1.0 - r2 * 0.6) * v_a * u_alpha;
writeDepth(v_depth);
o_col = outc(v_col * 0.5, a);
}`;
SH.VS_MARK = SH.HEAD + `
in vec3 a_pos; in vec4 a_col; uniform float u_dpr; uniform float u_size;
out vec4 v_col; out float v_depth;
void main(){ vec4 vp = u_view * vec4(a_pos, 1.0); v_depth = -vp.z; gl_Position = u_proj * vp; gl_PointSize = u_size * u_dpr; v_col = a_col; }`;
SH.FS_MARK = SH.HEAD + SH.FRAG + `
in vec4 v_col; in float v_depth;
void main(){
vec2 d = gl_PointCoord * 2.0 - 1.0; float r = length(d);
if (r > 1.0 || v_col.a < 0.01) discard;
float a = (1.0 - smoothstep(0.35, 0.6, r)) * v_col.a;
writeDepth(v_depth * 0.9995);
o_col = outc(v_col.rgb, a);
}`;
SH.VS_CRAFT = SH.HEAD + `
in vec3 a_pos; in vec3 a_nrm; in vec4 a_col; uniform mat4 u_model;
out vec4 v_col; out vec3 v_wn; out vec3 v_wp; out float v_depth;
void main(){ vec4 wp = u_model * vec4(a_pos, 1.0); v_wp = wp.xyz; v_wn = normalize(mat3(u_model) * a_nrm); v_col = a_col; vec4 vp = u_view * wp; v_depth = -vp.z; gl_Position = u_proj * vp; }`;
SH.FS_CRAFT = SH.HEAD + SH.FRAG + `
in vec4 v_col; in vec3 v_wn; in vec3 v_wp; in float v_depth;
uniform vec4 u_occ[4]; uniform int u_occN; uniform float u_sunR; uniform float u_ambient; uniform float u_hover;
uniform vec4 u_eng; uniform vec3 u_engCol; uniform float u_engR;
float eclipse(vec3 wp){
float f = 1.0; vec3 ts = u_sun - wp; float ds = length(ts); vec3 s = ts / ds; float as = u_sunR / ds;
for (int i = 0; i < 4; i++){
if (i >= u_occN) break;
vec3 c = u_occ[i].xyz - wp; float r = u_occ[i].w; float t = dot(c, s);
if (t <= 0.0 || t >= ds) continue;
float beta = length(c - s * t) / t, ao = r / t;
f *= 1.0 - (1.0 - smoothstep(abs(ao - as), ao + as, beta)) * clamp(sq(ao) / sq(as), 0.0, 1.0);
}
return f;
}
void main(){
vec3 N = normalize(v_wn); if (!gl_FrontFacing) N = -N;
vec3 L = normalize(u_sun - v_wp), Vv = normalize(-v_wp);
float ndl = dot(N, L);
float sh = eclipse(v_wp);
vec3 base = v_col.rgb;
vec3 H = normalize(L + Vv);
float sp = pow(max(dot(N, H), 0.0), 48.0) * v_col.a * step(0.0, ndl);
vec3 col = base * (max(ndl, 0.0) * 1.2 * sh + u_ambient + 0.03 + max(dot(N, Vv), 0.0) * 0.2) + vec3(1.0, 0.94, 0.82) * sp * sh * 1.2;
col += vec3(0.55, 0.7, 1.0) * pow(max(1.0 - max(dot(N, Vv), 0.0), 0.0), 3.0) * u_hover * 0.3;
if (u_eng.w > 0.0){ vec3 le = u_eng.xyz - v_wp; float de = length(le); col += base * u_engCol * u_eng.w * (0.35 + 0.65 * max(dot(N, le / de), 0.0)) / (1.0 + sq(de / u_engR)); }
writeDepth(v_depth);
o_col = outc(col, 1.0);
}`;
SH.VS_CRAFT_TEX = SH.VS_CRAFT.replace('in vec4 a_col;', 'in vec4 a_col; in vec2 a_uv;').replace('out vec4 v_col;', 'out vec4 v_col; out vec2 v_uv;').replace('v_col = a_col;', 'v_col = a_col; v_uv = a_uv;');
SH.FS_CRAFT_TEX = SH.FS_CRAFT.replace('in vec4 v_col;', 'in vec4 v_col; in vec2 v_uv; uniform sampler2D u_tex; uniform float u_hasTex; uniform sampler2D u_etex; uniform float u_eStr;')
.replace('  writeDepth(v_depth);\n  o_col = outc(col, 1.0);', '  if (u_eStr > 0.0) col += pow(texture(u_etex, v_uv).rgb, vec3(2.2)) * u_eStr;\n  writeDepth(v_depth);\n  o_col = outc(col, 1.0);')
.replace('vec3 base = v_col.rgb;', 'vec3 base = v_col.rgb * (u_hasTex > 0.5 ? pow(texture(u_tex, v_uv).rgb, vec3(2.2)) : vec3(1.0));');
SH.POST = `#version 300 es
precision highp float;
in vec2 v_uv; out vec4 o_col;
`;
SH.FS_DOWN = (thresh) => SH.POST + (thresh ? '#define THRESH\n' : '') + `
uniform sampler2D u_src; uniform vec2 u_texel; uniform float u_threshold;
vec3 S(vec2 o){
vec3 c = texture(u_src, v_uv + o * u_texel).rgb;
#ifdef THRESH
if (any(isnan(c)) || any(isinf(c))) return vec3(0.0);
c = min(max(c, vec3(0.0)), vec3(64.0));
#endif
return c;
}
void main(){
vec3 a = S(vec2(-2.0, 2.0)), b = S(vec2(0.0, 2.0)), c = S(vec2(2.0, 2.0));
vec3 d = S(vec2(-2.0, 0.0)), e = S(vec2(0.0)), f = S(vec2(2.0, 0.0));
vec3 g = S(vec2(-2.0, -2.0)), h = S(vec2(0.0, -2.0)), i = S(vec2(2.0, -2.0));
vec3 j = S(vec2(-1.0, 1.0)), k = S(vec2(1.0, 1.0)), l = S(vec2(-1.0, -1.0)), m = S(vec2(1.0, -1.0));
vec3 col = e * 0.125 + (a + c + g + i) * 0.03125 + (b + d + f + h) * 0.0625 + (j + k + l + m) * 0.125;
#ifdef THRESH
float br = max(col.r, max(col.g, col.b)); float knee = u_threshold * 0.5;
float soft = clamp(br - u_threshold + knee, 0.0, 2.0 * knee); soft = soft * soft / (4.0 * knee + 1e-5);
float w = max(soft, br - u_threshold) / max(br, 1e-5);
col = min(col * w, vec3(40.0));
#endif
o_col = vec4(col, 1.0);
}`;
SH.FS_UP = SH.POST + `
uniform sampler2D u_src; uniform vec2 u_texel;
void main(){
vec2 t = u_texel; vec3 c = vec3(0.0);
c += texture(u_src, v_uv + vec2(-t.x, t.y)).rgb + texture(u_src, v_uv + vec2(t.x, t.y)).rgb + texture(u_src, v_uv + vec2(-t.x, -t.y)).rgb + texture(u_src, v_uv + vec2(t.x, -t.y)).rgb;
c += 2.0 * (texture(u_src, v_uv + vec2(0.0, t.y)).rgb + texture(u_src, v_uv + vec2(0.0, -t.y)).rgb + texture(u_src, v_uv + vec2(t.x, 0.0)).rgb + texture(u_src, v_uv + vec2(-t.x, 0.0)).rgb);
c += 4.0 * texture(u_src, v_uv).rgb;
o_col = vec4(c / 16.0, 1.0);
}`;
SH.FS_COMP = SH.POST + `
uniform sampler2D u_scene; uniform sampler2D u_bloom; uniform float u_bloomK; uniform float u_exposure; uniform float u_inScale; uniform vec2 u_res; uniform float u_vignette;
vec3 aces(vec3 x){ const float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14; return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0); }
void main(){
vec3 c = texture(u_scene, v_uv).rgb * u_inScale;
vec3 bl = texture(u_bloom, v_uv).rgb * u_inScale;
if (any(isnan(c)) || any(isinf(c))) c = vec3(0.0);
if (any(isnan(bl)) || any(isinf(bl))) bl = vec3(0.0);
c = max(c, vec3(0.0)); bl = max(bl, vec3(0.0));
c = (c + bl * u_bloomK) * u_exposure;
c = aces(c);
vec2 q = v_uv - 0.5; c *= 1.0 - u_vignette * dot(q, q) * 1.3;
c = pow(c, vec3(1.0 / 2.2));
float nz = fract(sin(dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
c += (nz - 0.5) / 255.0;
o_col = vec4(c, 1.0);
}`;
const S = {
scale: 'visual', quality: 'alta', ambient: 0.02, rotCap: true, mode: 'explore',
layers: {
orbits: true, moonOrbits: true, cometOrbits: false, asteroids: true, kuiper: true, moons: true,
namesPlanets: true, namesMoons: true, namesSmall: true, constNames: false, constLines: false,
orbitInfo: false, deepSky: true, stars: true, markers: true, craft: true, namesCraft: true,
hz: false,          // capa científica: zona habitable
},
hzModel: 'conservador',
};
const COLORS = { orbitPlanet: '#8fb3e6', orbitDwarf: '#c2ae8c', orbitAst: '#a39886', orbitComet: '#7fd0ff', orbitMoon: '#8095b5', orbitCraft: '#86d3c6', select: '#f2c879' };
const lin = h => srgbToLin(hexToRgb(h));
const World = {
rb: [], byId: {}, sun: null, hd: {}, jd: 0, cam: [0, 0, 0], view: null, proj: null, meshes: {}, P: {},
sky: null, stars: null, constLines: null, belts: {}, markers: null, comets: [], landTex: null,
shadowSunR: 1, frame: 0,
async build(progress) {
const gl = GLX.gl;
const step = async (msg, f) => { progress(msg); await new Promise(r => requestAnimationFrame(() => r())); f(); };
await step('Compilando sombreadores', () => this.buildPrograms());
await step('Generando geometría', () => this.buildMeshes());
await step('Dibujando continentes', () => { this.landTex = GLX.texFromCanvas(EarthMap.build()); });
await step('Pintando la Vía Láctea', () => this.bakeSky());
await step('Situando estrellas y constelaciones', () => this.buildStars());
this.ptex = {}; this.texState = {};       // texturas y modelos detallados: bajo demanda (ver texTick y PackedModels.tick)
await step('Calculando órbitas', () => { this.buildBodies(); this.applyPlanetTextures(); this.buildSystems(); });
await step('Sembrando asteroides', () => this.buildBelts());
await step('Preparando cometas', () => this.buildComets());
gl.bindVertexArray(null);
},
buildPrograms() {
const P = this.P, p = (n, vs, fs) => (P[n] = GLX.program(n, vs, fs));
const feats = new Set();
BODIES.forEach(b => { if (b.vis.style === 'rock') feats.add(b.vis.feat); });
feats.forEach(f => p('rock_' + f, SH.VS_BODY, SH.FS_ROCK(f)));
['jupiter', 'saturno', 'urano', 'neptuno'].forEach(f => p('gas_' + f, SH.VS_BODY, SH.FS_GAS(f)));
p('earth', SH.VS_BODY, SH.FS_EARTH); p('clouds', SH.VS_BODY, SH.FS_CLOUDS); p('venus', SH.VS_BODY, SH.FS_VENUS);
p('atm', SH.VS_BODY, SH.FS_ATM); p('sun', SH.VS_BODY, SH.FS_SUN);
['saturn', 'uranus', 'faint', 'neptune'].forEach(t => p('ring_' + t, SH.VS_RING, SH.FS_RING(t)));
p('corona', SH.VS_BILL, SH.FS_CORONA); p('glow', SH.VS_BILL, SH.FS_GLOW); p('tail', SH.VS_TAIL, SH.FS_TAIL);
p('orbit', SH.VS_ORBIT, SH.FS_ORBIT); p('stars', SH.VS_STARS, SH.FS_STARS); p('linesky', SH.VS_LINESKY, SH.FS_LINESKY);
p('skybake', SH.VS_FULL, SH.FS_SKYBAKE); p('sky', SH.VS_FULL, SH.FS_SKY);
p('belt', SH.VS_BELT, SH.FS_BELT); p('mark', SH.VS_MARK, SH.FS_MARK); p('craft', SH.VS_CRAFT, SH.FS_CRAFT); p('hz', SH.VS_HZ, SH.FS_HZ); p('bh', SH.VS_FULL, SH.FS_BH); ['rock', 'luna', 'earth', 'venus', 'titan', 'gas'].forEach(f => p('ptex_' + f, SH.VS_BODY, SH.FS_PTEX(f))); p('ptex_gasring', SH.VS_BODY, SH.FS_PTEX('gas').replace('#define F_GAS', '#define F_GAS\n#define F_RINGTEX')); p('sunTex', SH.VS_BODY, SH.FS_SUNTEX); p('ring_saturnTex', SH.VS_RING, SH.FS_RING('saturnTex')); p('craftTex', SH.VS_CRAFT_TEX, SH.FS_CRAFT_TEX);
this.post = { progs: { downT: GLX.program('downT', SH.VS_FULL, SH.FS_DOWN(true)), down: GLX.program('down', SH.VS_FULL, SH.FS_DOWN(false)), up: GLX.program('up', SH.VS_FULL, SH.FS_UP), comp: GLX.program('comp', SH.VS_FULL, SH.FS_COMP) } };
},
buildMeshes() {
const gl = GLX.gl, M = this.meshes;
const sph = (a, b) => { const g = Geo.sphere(a, b); return GLX.mesh({ a_pos: { data: g.pos, size: 3 }, a_nrm: { data: g.nrm, size: 3 }, a_dir: { data: g.dir, size: 3 } }, g.idx); };
M.hi = sph(192, 96); M.mid = sph(96, 48); M.low = sph(36, 18);
const rg = (a, b) => { const g = Geo.ring(a, b, 256, 6); return GLX.mesh({ a_pos: { data: g.pos, size: 3 } }, g.idx); };
M.quad = GLX.mesh({ a_pos: { data: Geo.quad(), size: 3 } }, null, gl.TRIANGLE_STRIP);
M.rings = {};
BODIES.forEach(b => { if (b.vis.rings) M.rings[b.id] = rg(b.vis.rings.in, b.vis.rings.out); });
const n = 256, cp = new Float32Array((n + 1) * 3), cf = new Float32Array(n + 1);
for (let i = 0; i <= n; i++) { const a = TAU * i / n; cp[i * 3] = Math.cos(a); cp[i * 3 + 2] = -Math.sin(a); cf[i] = i / n; }
M.circle = GLX.mesh({ a_pos: { data: cp, size: 3 }, a_frac: { data: cf, size: 1 } }, null, gl.LINE_STRIP);
},
galBasis() {
const gnp = V.norm(Astro.radecToScene(192.85948, 27.12825));
let gc = Astro.radecToScene(266.405, -28.93617);
gc = V.norm(V.sub(gc, V.scale(gnp, V.dot(gc, gnp))));
return { gnp, gc, gy: V.cross(gnp, gc) };
},
bakeSky() {
const gl = GLX.gl, size = Math.min({ basic: 512, high: 1024, max: 2048 }[Settings.state.graphics.sky] || 1024, gl.getParameter(gl.MAX_CUBE_MAP_TEXTURE_SIZE));
const tex = gl.createTexture(); this.skySize = size;
gl.bindTexture(gl.TEXTURE_CUBE_MAP, tex);
for (let i = 0; i < 6; i++) {
if (GLX.hdr) gl.texImage2D(gl.TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, gl.RGBA16F, size, size, 0, gl.RGBA, gl.HALF_FLOAT, null);
else gl.texImage2D(gl.TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, gl.RGBA8, size, size, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
}
gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
const fbo = gl.createFramebuffer(), pr = this.P.skybake, G = this.galBasis();
const pDir = [], pCol = [], pAx = [];
DEEP_SKY.forEach(d => {
const ra = d[1] * 15 * DEG, dec = d[2] * DEG, pa = d[7] * DEG;
const dir = Astro.radecToScene(d[1] * 15, d[2]);
const e = [-Math.sin(ra), Math.cos(ra), 0], nn = [-Math.sin(dec) * Math.cos(ra), -Math.sin(dec) * Math.sin(ra), Math.cos(dec)];
const ax = V.norm(Astro.eclToScene(Astro.eqToEcl(V.add(V.scale(nn, Math.cos(pa)), V.scale(e, Math.sin(pa))))));
pDir.push(dir[0], dir[1], dir[2], d[3] * DEG); pCol.push(d[4][0], d[4][1], d[4][2], d[5]); pAx.push(ax[0], ax[1], ax[2], d[6]);
});
while (pDir.length < 48) { pDir.push(0, 0, 0, 0); pCol.push(0, 0, 0, 0); pAx.push(1, 0, 0, 1); }
GLX.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'none' });
GLX.use(pr);
GLX.setAll(pr, { u_size: size, u_gnp: G.gnp, u_gc: G.gc, u_gy: G.gy, u_pDir: new Float32Array(pDir), u_pCol: new Float32Array(pCol), u_pAx: new Float32Array(pAx), u_pN: DEEP_SKY.length });
gl.bindVertexArray(GLX.fsTri);
gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
gl.viewport(0, 0, size, size);
for (let i = 0; i < 6; i++) {
gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_CUBE_MAP_POSITIVE_X + i, tex, 0);
GLX.set(pr, 'u_face', i);
gl.drawArrays(gl.TRIANGLES, 0, 3);
}
gl.bindFramebuffer(gl.FRAMEBUFFER, null);
gl.deleteFramebuffer(fbo);
if (this.sky) gl.deleteTexture(this.sky);
this.sky = tex;
},
buildStars() {
const G = this.galBasis(), rng = makeRng(42), pos = [], col = [];
const tcol = t => { // color aproximado de cuerpo negro (sRGB → lineal)
const k = t / 100; let r, g, b;
r = k <= 66 ? 255 : 329.7 * Math.pow(k - 60, -0.1332);
g = k <= 66 ? 99.47 * Math.log(k) - 161.1 : 288.1 * Math.pow(k - 60, -0.0755);
b = k >= 66 ? 255 : k <= 19 ? 0 : 138.5 * Math.log(k - 10) - 305;
return srgbToLin([clamp(r, 0, 255) / 255, clamp(g, 0, 255) / 255, clamp(b, 0, 255) / 255]);
};
const N = 11000; this.starN = N;            // la cantidad visible la decide el ajuste de partículas
for (let i = 0; i < N; i++) {
let b;
if (rng() < 0.42) b = gauss(rng) * 13 * DEG; else b = Math.asin(rng() * 2 - 1);
const l = rng() * TAU;
const d = V.add(V.add(V.scale(G.gc, Math.cos(b) * Math.cos(l)), V.scale(G.gy, Math.cos(b) * Math.sin(l))), V.scale(G.gnp, Math.sin(b)));
const mag = 7.2 - 5.0 * Math.pow(rng(), 2.6);
const r = rng(), T = r < 0.12 ? 3200 + rng() * 800 : r < 0.45 ? 4200 + rng() * 1600 : r < 0.85 ? 5800 + rng() * 2400 : 8500 + rng() * 16000;
const c = tcol(T);
pos.push(d[0], d[1], d[2]); col.push(c[0], c[1], c[2], mag);
}
this.namedStars = [];
for (const k in STARS) {
const s = STARS[k], d = Astro.radecToScene(s[1] * 15, s[2]);
const ci = { betelgeuse: 3300, antares: 3400, aldebaran: 3900, arcturus: 4300, pollux: 4700, rigel: 12000, spica: 22000, vega: 9600, sirius: 9900 }[k] || 6500 + ((k.length * 977) % 5000);
const c = tcol(ci);
pos.push(d[0], d[1], d[2]); col.push(c[0], c[1], c[2], s[3]);
this.namedStars.push({ id: k, name: s[0], dir: d });
}
this.stars = GLX.mesh({ a_pos: { data: new Float32Array(pos), size: 3 }, a_col: { data: new Float32Array(col), size: 4 } }, null, GLX.gl.POINTS);
const lp = [];
this.constellations = CONSTELLATIONS.map(c => {
let cen = [0, 0, 0];
const ids = new Set();
c.lines.forEach(([a, b]) => {
const A = Astro.radecToScene(STARS[a][1] * 15, STARS[a][2]), B = Astro.radecToScene(STARS[b][1] * 15, STARS[b][2]);
lp.push(...A, ...B); ids.add(a); ids.add(b);
});
ids.forEach(id => { cen = V.add(cen, Astro.radecToScene(STARS[id][1] * 15, STARS[id][2])); });
return { id: c.id, name: c.name, dir: V.norm(cen) };
});
this.constLines = GLX.mesh({ a_pos: { data: new Float32Array(lp), size: 3 } }, null, GLX.gl.LINES);
},
makeBody(def) {
const jd0 = Astro.jdFromDate(new Date());
const rb = { def, id: def.id, isSun: def.type === 'star', R: def.R, children: [], hover: 0, rotPhase: null,
helio: [0, 0, 0], posS: [0, 0, 0], rS: 1, basis: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], model: new Float64Array(16),
proj: { x: 0, y: 0, w: 0, rpx: 0, on: false }, color: lin(def.color || '#cccccc'), occ: [], occBuf: new Float32Array(16) };
rb.frame = def.pole ? Astro.poleFrame(def.pole[0], def.pole[1]) : { P: [0, 1, 0], Q: [1, 0, 0], PxQ: V.cross([0, 1, 0], [1, 0, 0]) };
rb.rate = def.rotH ? 360 * 24 / def.rotH : 0;
const o = def.orbit;
if (o && o.t === 'kep') {
if (o.tp) o.tpJD = Astro.jdFromISO(o.tp);
else if (o.fit) this.fitOrbit(o);
}
const v = def.vis;
if (def.shape || v.irregular || v.feat === 'bilobe') {
const g = Geo.irregular(def.shape, v.irregular || 0.05, v.seed || 1, v.feat === 'bilobe');
rb.mesh = GLX.mesh({ a_pos: { data: g.pos, size: 3 }, a_nrm: { data: g.nrm, size: 3 }, a_dir: { data: g.dir, size: 3 } }, g.idx);
rb.irregular = true;
}
if (v.style === 'craft') {
rb.mesh = CraftModels.mesh(def.model); rb.isCraft = true;
if (def.model === 'jwst' && this.jwstHD) rb.lodMesh = CraftModels.mesh('jwst-lite');
rb.vFrom = def.validFrom ? Astro.jdFromISO(def.validFrom) : null; rb.vTo = def.validTo ? Astro.jdFromISO(def.validTo) : null;
if (o && o.t === 'drift') { o.dir = V.norm(Astro.eqToEcl(Astro.radecVec(o.ra, o.dec))); o.t0JD = Astro.jdFromISO(o.t0); }
}
rb.prog = v.style === 'rock' ? this.P['rock_' + v.feat] : v.style === 'gas' ? this.P['gas_' + v.feat] : this.P[v.style];
rb.u = {};
if (v.c1) { rb.u.u_c1 = lin(v.c1); rb.u.u_c2 = lin(v.c2); rb.u.u_c3 = lin(v.c3); rb.u.u_crater = v.crater; rb.u.u_seed = v.seed; }
rb.u.u_bump = v.style === 'earth' ? 0.012 : v.style === 'rock' ? (rb.irregular ? 0.05 : 0.03) : 0;
rb.u.u_atm = v.atm ? [...v.atm.color.map(x => x * 0.9), v.atm.k] : [0, 0, 0, 0];
if (v.style === 'earth') rb.u.u_land = this.landTex;
return rb;
},
buildBodies() {
const jd0 = Astro.jdFromDate(new Date());
BODIES.forEach(def => { const rb = this.makeBody(def); this.rb.push(rb); this.byId[rb.id] = rb; });
this.rb.forEach(rb => { if (rb.def.parent) { rb.parent = this.byId[rb.def.parent]; rb.parent.children.push(rb); } });
this.sun = this.byId.sol;
this.rb.forEach(rb => {
const o = rb.def.orbit;
if (o && (o.t === 'moon' || o.t === 'leo')) rb.oframe = Astro.orbitFrame(rb.parent.frame.P, rb.parent.frame.Q, o.i, o.node);
});
this.rb.forEach(rb => {
if (rb.isSun) return;
if (rb.isCraft) rb.occ = rb.parent && !rb.parent.isSun ? [rb.parent] : [];
else if (rb.def.type === 'moon') rb.occ = [rb.parent].concat(rb.parent.children.filter(c => c !== rb && c.def.type === 'moon').sort((a, b) => b.R - a.R).slice(0, 3));
else rb.occ = rb.children.filter(c => c.def.type === 'moon').sort((a, b) => b.R - a.R).slice(0, 4);
});
this.rb.forEach(rb => {
const o = rb.def.orbit; if (!o || (o.t !== 'jpl' && o.t !== 'kep')) return;
const n = rb.def.type === 'comet' ? 1024 : (o.e > 0.3 ? 720 : 512);
const s = Astro.sampleOrbit(o, jd0, n);
rb.orbitPts = s.pts.map(p => Astro.eclToScene(p));
rb.orbitBuf = new Float32Array(rb.orbitPts.length * 3);
rb.orbitMesh = GLX.mesh({ a_pos: { data: rb.orbitBuf, size: 3 }, a_frac: { data: new Float32Array(s.fr), size: 1 } }, null, GLX.gl.LINE_STRIP, true);
});
const mk = this.rb.length;
this.markBuf = { pos: new Float32Array(mk * 3), col: new Float32Array(mk * 4) };
this.markers = GLX.mesh({ a_pos: { data: this.markBuf.pos, size: 3 }, a_col: { data: this.markBuf.col, size: 4 } }, null, GLX.gl.POINTS, true);
this.remapOrbits();
},
drawGroups(mesh, pr) {
for (const g of mesh.groups) {
pr.tex = 0;
GLX.set(pr, 'u_tex', g.tex || GLX.whiteTex()); GLX.set(pr, 'u_hasTex', g.tex ? 1 : 0);
GLX.set(pr, 'u_etex', g.etex || GLX.whiteTex()); GLX.set(pr, 'u_eStr', g.etex ? g.estr : 0);
GLX.drawRange(mesh, pr, g.first, g.count);
}
},
requestTex(id, visible) {
if (typeof PLANET_TEX === 'undefined' || !PLANET_TEX[id] || this.texState[id]) return;
this.texState[id] = 'loading'; Assets.changed();
if (visible) Assets.busy('t:' + id, 'Preparando objeto…');
if (ASSET_MANIFEST.bodyMeshes && ASSET_MANIFEST.bodyMeshes[id]) PackedModels.loadBody(id);   // forma real del cuerpo
(async () => {
const e = PLANET_TEX[id], map = await Assets.image(e.map, e.mapAlt);
if (!map) throw new Error('sin formato compatible');
const t = { map: await GLX.texFromURI(map, true) };
if (e.nmap) t.nmap = await GLX.texFromURI(e.nmap, true);
if (e.ring) t.ring = await GLX.texFromURI(e.ring.uri, true);
this.ptex[id] = t; this.applyPlanetTextures(id); this.texState[id] = 'ok'; Assets.changed();
if (UI.sel === this.byId[id]) UI.renderInfo(UI.sel);          // la ficha abierta muestra ya el crédito del modelo
})().catch(err => { this.texState[id] = 'fail'; Assets.changed(); if (!/formato/.test(err.message)) console.warn('Textura no disponible: ' + id, err); })
.finally(() => Assets.done('t:' + id));
},
texTick() {
if (typeof PLANET_TEX === 'undefined' || !this.texState) return;
let pending = 0;
for (const id in PLANET_TEX) {
const rb = this.byId[id]; if (!rb || this.texState[id]) { if (this.texState[id] === 'loading') pending++; continue; }
const focus = UI.sel === rb || UI.tourRB === rb || (Flight.on && Flight.target === rb);
if (rb.isSun || id === 'tierra' || focus || (rb.drawOn && rb.proj.rpx > 2.5)) { this.requestTex(id, focus); pending++; }
}
this.idleT = (this.idleT || 0) + 1;
if (!pending && this.idleT > 24) {
const order = ['marte', 'jupiter', 'saturno', 'venus', 'mercurio', 'luna', 'urano', 'neptuno', 'pluton', 'ganimedes', 'titan'];
const next = order.find(id => PLANET_TEX[id] && !this.texState[id]); if (next) this.requestTex(next, false);
}
},
applyPlanetTextures(only) {
const feat = { mercurio: 'rock', marte: 'rock', pluton: 'rock', luna: 'luna', ganimedes: 'luna', titan: 'titan', tierra: 'earth', venus: 'venus', jupiter: 'gas', saturno: 'gasring', urano: 'gas', neptuno: 'gas',
io: 'rock', calisto: 'luna', rea: 'luna', oberon: 'luna', ariel: 'luna', miranda: 'luna', fobos: 'rock', deimos: 'rock',
ceres: 'rock', haumea: 'luna', makemake: 'luna', sedna: 'rock', quaoar: 'luna' };
for (const id of only ? [only] : Object.keys(this.ptex)) {
const rb = this.byId[id], t = this.ptex[id]; if (!rb || !t) continue;
this.hd[id] = true;
if (rb.isSun) { this.sunTex = t.map; continue; }
rb.prog = this.P['ptex_' + feat[id]];
rb.u.u_map = t.map; rb.u.u_lunar = { mercurio: 0.75, luna: 0.75, marte: 0.3, pluton: 0.4, ganimedes: 0.6, io: 0.4, calisto: 0.6, rea: 0.7, oberon: 0.6, ariel: 0.6, miranda: 0.6, fobos: 0.75, deimos: 0.75, ceres: 0.7, haumea: 0.5, makemake: 0.5, sedna: 0.4, quaoar: 0.5 }[id] || 0;
rb.u.u_bump = feat[id] === 'rock' ? 0.012 : 0;
if (t.nmap) { rb.u.u_nmap = t.nmap; rb.u.u_nScale = PLANET_TEX[id].nScale; }
else { rb.u.u_nmap = t.map; rb.u.u_nScale = 0; }   // sin mapa de relieve propio: relieve nulo (no hereda el del último cuerpo dibujado)
if (id === 'tierra') { rb.u.u_land = this.landTex; this.earthTex = true; }
if (t.ring) { this.ringTex = t.ring; rb.u.u_ringTex = t.ring; rb.u.u_rr2 = [PLANET_TEX[id].ring.rin, PLANET_TEX[id].ring.rout]; }
}
},
fitOrbit(o) {
const jd = Astro.jdFromISO(o.fit.date), n = GAUSS_K_DEG / Math.pow(o.a, 1.5);
let target;
if (o.fit.earth) target = Astro.helioPos(BODY.tierra.orbit, jd);
else target = V.norm(Astro.eqToEcl(Astro.radecVec(o.fit.ra, o.fit.dec)));
const E = o.fit.earth ? null : Astro.helioPos(BODY.tierra.orbit, jd);
const score = M => {
const k = Astro.keplerE(M * DEG, o.e), p = Astro.elemToEcl(o.a, o.e, o.i, o.node, o.peri, k);
if (o.fit.earth) return V.dist(p, target);
return -V.dot(V.norm(V.sub(p, E)), target);
};
let best = 0, bs = Infinity;
for (let M = 0; M < 360; M += 0.5) { const s = score(M); if (s < bs) { bs = s; best = M; } }
for (let st = 0.25; st > 1e-4; st *= 0.5) for (const d of [-st, st]) { const s = score(best + d); if (s < bs) { bs = s; best += d; } }
o.tpJD = jd - best / n;
},
remapOrbits() {
this.rb.forEach(rb => {
if (!rb.orbitMesh) return;
const b = rb.orbitBuf;
rb.orbitPts.forEach((p, i) => { const m = ScaleState.mapVec(p); b[i * 3] = m[0]; b[i * 3 + 1] = m[1]; b[i * 3 + 2] = m[2]; });
GLX.update(rb.orbitMesh, 'a_pos', b);
});
},
buildBelts() {
const rng = makeRng(2024), e1 = [], e2 = [], k1 = [], k2 = [];
const push = (A, B, a, e, i, node, peri, M0, size, tint) => { A.push(a, e, i * DEG, node * DEG); B.push(peri * DEG, M0 * DEG, size, tint); };
const gaps = [[2.502, 0.03], [2.825, 0.02], [2.958, 0.014], [3.279, 0.04]];
const nMain = 8000;
for (let c = 0; c < nMain;) {
const a = 2.06 + 1.26 * Math.pow(rng(), 0.9);
if (gaps.some(([g, w]) => Math.abs(a - g) < w * (0.6 + rng() * 0.6))) continue;
if (a > 3.0 && rng() < 0.45) continue;
push(e1, e2, a, Math.min(0.32, Math.abs(gauss(rng)) * 0.09 + 0.02), Math.abs(gauss(rng)) * 7.5, rng() * 360, rng() * 360, rng() * 360, 0.6 + rng() * 0.8, rng() * 0.5);
c++;
}
const LJ = 34.40;
for (let c = 0; c < 1500; c++) {
const node = rng() * 360, peri = rng() * 360, lam = LJ + (rng() < 0.6 ? 60 : -60) + gauss(rng) * 11;
push(e1, e2, 5.2026 + gauss(rng) * 0.04, Math.abs(gauss(rng)) * 0.06, Math.abs(gauss(rng)) * 12, node, peri, lam - node - peri, 0.7 + rng() * 0.6, 0.35 + rng() * 0.15);
}
for (let c = 0; c < 260; c++) push(e1, e2, 0.95 + rng() * 1.6, 0.1 + rng() * 0.5, Math.abs(gauss(rng)) * 10, rng() * 360, rng() * 360, rng() * 360, 0.5 + rng() * 0.4, 0.2);
const nK = 5500;
for (let c = 0; c < nK; c++) {
const r = rng();
if (r < 0.55) push(k1, k2, 42.2 + rng() * 5.5, rng() * 0.09, Math.abs(gauss(rng)) * 3.5, rng() * 360, rng() * 360, rng() * 360, 1 + rng(), 0.55 + rng() * 0.45);
else if (r < 0.8) push(k1, k2, 39.45 + gauss(rng) * 0.15, 0.1 + rng() * 0.2, Math.abs(gauss(rng)) * 10, rng() * 360, rng() * 360, rng() * 360, 1 + rng(), 0.55 + rng() * 0.45);
else push(k1, k2, 32 + rng() * 60, 0.2 + rng() * 0.4, Math.abs(gauss(rng)) * 15, rng() * 360, rng() * 360, rng() * 360, 1 + rng(), 0.6 + rng() * 0.4);
}
const shuffle = (A, B) => { const n = A.length / 4, idx = [...Array(n).keys()]; for (let i = n - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [idx[i], idx[j]] = [idx[j], idx[i]]; }
const A2 = new Float32Array(A.length), B2 = new Float32Array(B.length); idx.forEach((o, i) => { for (let k = 0; k < 4; k++) { A2[i * 4 + k] = A[o * 4 + k]; B2[i * 4 + k] = B[o * 4 + k]; } }); return [A2, B2]; };
const mk = (A, B) => { const [A2, B2] = shuffle(A, B); return GLX.mesh({ a_e1: { data: A2, size: 4 }, a_e2: { data: B2, size: 4 } }, null, GLX.gl.POINTS); };
this.belts.ast = mk(e1, e2); this.belts.kui = mk(k1, k2);
},
buildComets() {
const gl = GLX.gl, N = 48;
this.rb.filter(r => r.def.type === 'comet').forEach((rb, i) => {
const mk = () => GLX.mesh({ a_pos: { data: new Float32Array((N + 1) * 2 * 3), size: 3 }, a_uv: { data: new Float32Array((N + 1) * 2 * 2), size: 2 } }, null, gl.TRIANGLE_STRIP, true);
const c = { rb, ion: mk(), dust: mk(), N, act: 0, seed: i * 3.7 };
rb.comet = c; this.comets.push(c);
});
},
update(jd, dtReal, speed) {
if (this.system !== 'solar') return this.updateSystem(jd, dtReal, speed);   // otros sistemas: su propia actualización
this.jd = jd;
const d = jd - J2000;
const sc = ScaleState;
if (sc.t < 1) { sc.t = Math.min(1, sc.t + dtReal / 1.6); this.remapOrbits(); }
for (const rb of this.rb) {
const o = rb.def.orbit;
rb.rS = sc.body(rb.R, rb.isSun);
if (rb.isSun) { rb.helio = [0, 0, 0]; rb.posS = [0, 0, 0]; continue; }
if (o.t === 'jpl' || o.t === 'kep') {
rb.helio = Astro.helioPos(o, jd);
rb.posS = sc.mapVec(Astro.eclToScene(rb.helio));
rb.orbitFrac = Astro.elementsAt(o, jd).M / 360;
}
}
for (const rb of this.rb) {
const o = rb.def.orbit; if (!o || (o.t !== 'moon' && o.t !== 'luna')) continue;
const par = rb.parent; let dirS, rKm;
if (o.t === 'luna') {
const L = Astro.lunaGeo(jd);
const v = Astro.eclToScene(L.pos); rKm = V.len(v); dirS = V.scale(v, 1 / rKm);
rb.oframe = Astro.orbitFrame([0, 1, 0], [1, 0, 0], 5.145, L.node);
rb.orbitFrac = L.frac;
} else {
const th = (o.L0 + 360 * d / o.P) * DEG;
dirS = V.norm(V.add(V.scale(rb.oframe.u, Math.cos(th)), V.scale(rb.oframe.v, Math.sin(th))));
rKm = o.a; rb.orbitFrac = ((th / TAU) % 1 + 1) % 1;
}
rb.moonDist = sc.moon(rKm, par.R, sc.bodyA(par.R), sc.bodyB(par.R));
rb.moonDistMean = sc.moon(o.a, par.R, sc.bodyA(par.R), sc.bodyB(par.R));
rb.posS = V.add(par.posS, V.scale(dirS, rb.moonDist));
rb.helio = V.add(par.helio, V.scale(Astro.sceneToEcl(dirS), rKm / AU_KM));
rb.moonDir = dirS;
}
for (const rb of this.rb) {
if (!rb.isCraft) continue;
const o = rb.def.orbit;
rb.hidden = !!(rb.vFrom && (jd < rb.vFrom || (rb.vTo && jd > rb.vTo)));
if (o.t === 'leo') {
const par = rb.parent, th = (o.L0 + 360 * d / o.P) * DEG;
const dirS = V.norm(V.add(V.scale(rb.oframe.u, Math.cos(th)), V.scale(rb.oframe.v, Math.sin(th))));
const map = k => k === 'real' ? SCALES.real.moon(o.a) : SCALES[k].body(par.R) * (1 + (o.a / par.R - 1) * 2.5);
rb.moonDist = rb.moonDistMean = sc.blend(map(sc.from), map(sc.to));
rb.posS = V.add(par.posS, V.scale(dirS, rb.moonDist));
rb.helio = V.add(par.helio, V.scale(Astro.sceneToEcl(dirS), o.a / AU_KM));
rb.moonDir = dirS; rb.orbitFrac = ((th / TAU) % 1 + 1) % 1;
} else if (o.t === 'lpoint') {
const E = rb.parent, s = V.norm(E.helio), z = [0, 0, 1], w = V.norm(V.cross(s, z));
const base = V.scale(s, (o.L === 2 ? 1 : -1) * o.dist);
const th = (o.phase + 360 * d / o.P) * DEG;
const off = V.add(base, V.add(V.scale(z, o.halo[0] * Math.sin(th)), V.scale(w, o.halo[1] * Math.cos(th))));
const dk = V.len(off), map = km => sc.moon(km, E.R, sc.bodyA(E.R), sc.bodyB(E.R));
rb.helio = V.add(E.helio, V.scale(off, 1 / AU_KM));
rb.posS = V.add(E.posS, V.scale(Astro.eclToScene(V.scale(off, 1 / dk)), map(dk)));
const k = map(o.dist) / o.dist;
rb.halo = { c: V.scale(Astro.eclToScene(V.norm(base)), map(o.dist)), X: V.scale(Astro.eclToScene(w), o.halo[1] * k), Z: V.scale(Astro.eclToScene(z), -o.halo[0] * k), N: Astro.eclToScene(s) };
rb.orbitFrac = ((th / TAU) % 1 + 1) % 1;
} else if (o.t === 'drift') {
rb.helio = V.scale(o.dir, o.r0 + o.rate * (jd - o.t0JD) / 365.25);
rb.posS = sc.mapVec(Astro.eclToScene(rb.helio));
}
let Y, X;
if (o.t === 'leo') { Y = rb.moonDir; X = V.norm(V.cross(rb.oframe.n, Y)); }
else {
Y = V.norm(Astro.eclToScene(V.scale(rb.helio, -1)));
if (rb.def.model === 'jwst') Y = V.scale(Y, -1);
X = V.norm(V.cross(Y, Math.abs(Y[1]) > 0.95 ? [1, 0, 0] : [0, 1, 0]));
}
rb.basis = [X, Y, V.cross(X, Y)];
}
for (const rb of this.rb) {
if (rb.isCraft) continue;
if (rb.def.locked) {
const Y = rb.oframe.n, X = V.norm(V.scale(rb.moonDir, -1)), Z = V.cross(X, Y);
rb.basis = [X, Y, Z];
} else {
const exact = (rb.def.W0 || 0) + rb.rate * d;
const visRate = Math.abs(rb.rate * speed) / 86400; // °/s reales
if (S.rotCap && visRate > 150 && speed !== 0) {
if (rb.rotPhase == null) rb.rotPhase = exact;
rb.rotPhase += Math.sign(rb.rate * speed) * 150 * dtReal;
rb.capped = true;
} else { rb.rotPhase = exact; rb.capped = false; }
rb.basis = Astro.basisFromFrame(rb.frame, rb.rotPhase);
}
}
for (const c of this.comets) {
const r = V.len(c.rb.helio);
c.act = clamp((4.8 - r) / 4.2, 0, 1);
c.act = Math.pow(c.act, 1.4) * Math.min(1.8, 1 / Math.max(r, 0.35));
c.r = r;
}
},
project(p, out) {
const v = this.view, P = this.proj, W = GLX.W / GLX.dpr, H = GLX.H / GLX.dpr;
const x = p[0] - this.cam[0], y = p[1] - this.cam[1], z = p[2] - this.cam[2];
const vx = v[0] * x + v[4] * y + v[8] * z, vy = v[1] * x + v[5] * y + v[9] * z, vz = v[2] * x + v[6] * y + v[10] * z;
const w = -vz;
out.w = w;
if (w <= 1e-12) { out.on = false; return out; }
out.x = ((P[0] * vx / w - P[8]) * 0.5 + 0.5) * W;
out.y = (1 - ((P[5] * vy / w - P[9]) * 0.5 + 0.5)) * H;
out.on = true;
return out;
},
pxPerUnit(w) { return this.proj[5] * (GLX.H / GLX.dpr) * 0.5 / w; },
render(t) {
const gl = GLX.gl, P = this.P, cam = this.cam, L = S.layers;
const nowc = performance.now(), dtc = this._lt ? Math.min(0.1, (nowc - this._lt) / 1000) : 0.016; this._lt = nowc;
this.cineK = (this.cineK || 0) + (((S.mode === 'cine' || UI.tourCine || (typeof Missions !== 'undefined' && Missions.run && Missions.run.cine)) ? 1 : 0) - (this.cineK || 0)) * Math.min(1, dtc * 2.2);
const vis1 = 1 - this.cineK;
const SOLAR = this.system === 'solar';
const sunRS = this.sun.rS;
const sunR = V.sub(this.sun.posS, cam);
this.shadowSunR = sunRS * ScaleState.get('shadowSun');
const glob = { u_view: this.view, u_proj: this.proj, u_depthK: 1e6, u_depthInv: 1 / Math.log2(1 + 1e10 * 1e6), u_time: t, u_sun: sunR, u_outScale: GLX.hdr ? 1 : 0.25 };
const useP = pr => { GLX.use(pr); GLX.setAll(pr, glob); return pr; };
GLX.beginScene();
if (L.deepSky && Settings.state.graphics.nebulae) {
GLX.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'none' });
const pr = useP(P.sky), v = this.view, f = 1 / this.proj[5];
GLX.setAll(pr, { u_sky: this.sky, u_camR: [v[0], v[4], v[8]], u_camU: [v[1], v[5], v[9]], u_camF: [-v[2], -v[6], -v[10]], u_tan: [1 / this.proj[0], 1 / this.proj[5]], u_off: [this.proj[8], this.proj[9]], u_shift: [this.proj[8], this.proj[9]], u_skyB: 1.0 });
gl.bindVertexArray(GLX.fsTri); gl.drawArrays(gl.TRIANGLES, 0, 3);
}
if (L.stars) {
GLX.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'add' });
const pr = useP(P.stars); GLX.setAll(pr, { u_dpr: GLX.dpr, u_bright: 0.55 * this.starDim });
const fr = Gfx.particles(), total = this.stars.count, nr = this.starN || total;
GLX.bindMesh(this.stars, pr); gl.drawArrays(gl.POINTS, 0, Math.round(nr * fr)); if (total > nr) gl.drawArrays(gl.POINTS, nr, total - nr);
}
if (L.constLines && vis1 > 0.003) {
GLX.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'alpha' });
const pr = useP(P.linesky); GLX.setAll(pr, { u_color: lin('#7f9cc8'), u_alpha: 0.32 * vis1 });
GLX.draw(this.constLines, pr);
}
if (!SOLAR && this.sun.isBH) this.drawBH(t, useP);
GLX.state({ depthTest: true, depthWrite: true, cull: 'back', blend: 'none' });
const vis = this.rb.filter(rb => rb.drawOn);
for (const rb of vis) {
if (rb.isCraft || rb.isBH || rb.occBH) continue;
const rel = V.sub(rb.posS, cam);
M4.fromBasis(rb.basis[0], rb.basis[1], rb.basis[2], rb.rS, rel, rb.model);
const pr = useP(rb.isSun ? (this.sunTex ? P.sunTex : P.sun) : rb.prog);
if (rb.isSun && this.sunTex) GLX.set(pr, 'u_map', this.sunTex);
let mesh = rb.mesh;
if (!mesh) { const L2 = Gfx.lod(); mesh = rb.proj.rpx > L2[0] ? this.meshes.hi : rb.proj.rpx > L2[1] ? this.meshes.mid : this.meshes.low; }
const ob = rb.occBuf; let n = 0; const occMax = Gfx.occluders();
for (const o of rb.occ) { if (n >= occMax) break; const r2 = V.sub(o.posS, cam); ob[n * 4] = r2[0]; ob[n * 4 + 1] = r2[1]; ob[n * 4 + 2] = r2[2]; ob[n * 4 + 3] = o.rS * (o.irregular ? 0.8 : 1); n++; }
GLX.setAll(pr, rb.u);
GLX.setAll(pr, { u_model: rb.model, u_occ: ob, u_occN: n, u_sunR: this.shadowSunR, u_ambient: S.ambient, u_hover: rb.hover, u_radius: rb.rS, u_detail: rb.proj.rpx > 90 ? 1 : 0, u_center: rel, u_sunI: lerp(4.2, 0.72, smoothstep(12, 260, this.sun.proj.rpx)) });
if (rb.id === 'saturno') GLX.setAll(pr, { u_ringN: rb.basis[1], u_ringOn: Gfx.ringShadow() ? 1 : 0 });
if (rb.u.u_nmap) { const B = rb.basis; GLX.set(pr, 'u_rot', [B[0][0], B[0][1], B[0][2], B[1][0], B[1][1], B[1][2], B[2][0], B[2][1], B[2][2]]); }
GLX.draw(mesh, pr);
rb.rel = rel;
}
GLX.state({ depthTest: true, depthWrite: true, cull: 'none', blend: 'none' });
for (const rb of vis) {
if (!rb.isCraft) continue;
const rel = V.sub(rb.posS, cam); rb.rel = rel;
M4.fromBasis(rb.basis[0], rb.basis[1], rb.basis[2], rb.rS, rel, rb.model);
const mesh = rb.lodMesh && rb.proj.rpx < Gfx.lod()[2] ? rb.lodMesh : rb.mesh;
const pr = useP(mesh.groups ? P.craftTex : P.craft), ob = rb.occBuf; let n = 0; const occMaxC = Gfx.occluders();
for (const o of rb.occ) { if (n >= occMaxC) break; const r2 = V.sub(o.posS, cam); ob[n * 4] = r2[0]; ob[n * 4 + 1] = r2[1]; ob[n * 4 + 2] = r2[2]; ob[n * 4 + 3] = o.rS; n++; }
GLX.setAll(pr, { u_model: rb.model, u_occ: ob, u_occN: n, u_sunR: this.shadowSunR, u_ambient: S.ambient, u_hover: rb.hover, u_eng: [0, 0, 0, 0] });
if (mesh.groups) this.drawGroups(mesh, pr); else GLX.draw(mesh, pr);
}
if (Flight.on) Flight.drawShip(useP);
if (!SOLAR && typeof Missions !== 'undefined' && Missions.run) Missions.draw(useP, t);
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'premul' });
for (const rb of vis) {
const R = rb.def.vis.rings; if (!R || rb.proj.rpx < 1.2) continue;
const tex = R.type === 'saturn' && this.ringTex;
const pr = useP(P['ring_' + R.type + (tex ? 'Tex' : '')]);
if (tex) GLX.setAll(pr, { u_ringTex: this.ringTex, u_rr2: rb.u.u_rr2 });
GLX.setAll(pr, { u_model: rb.model, u_center: rb.rel, u_radius: rb.rS, u_ringN: rb.basis[1], u_ambient: S.ambient * 0.5, u_alpha: R.alpha, u_rr: [R.in, R.out] });
GLX.draw(this.meshes.rings[rb.id], pr);
}
if (SOLAR) {
const earth = this.byId.tierra;
if (earth.drawOn && earth.proj.rpx > 3 && !this.earthTex) {
GLX.state({ depthTest: true, depthWrite: false, cull: 'back', blend: 'premul' });
const pr = useP(P.clouds), m = M4.fromBasis(earth.basis[0], earth.basis[1], earth.basis[2], earth.rS * 1.008, earth.rel);
const ob = earth.occBuf;
GLX.setAll(pr, { u_model: m, u_cloudT: t * 0.004, u_occ: ob, u_occN: earth.occ.length ? 1 : 0, u_sunR: this.shadowSunR, u_ambient: S.ambient * 0.5 });
GLX.draw(earth.proj.rpx > 160 ? this.meshes.hi : this.meshes.mid, pr);
}
}
GLX.state({ depthTest: true, depthWrite: false, cull: 'front', blend: 'add' });
if (Settings.state.graphics.atmospheres) for (const rb of vis) {
const A = rb.def.vis.atm; if (!A || rb.proj.rpx < 2) continue;
const H = rb.rS * A.h, pr = useP(P.atm);
const m = M4.fromBasis([1, 0, 0], [0, 1, 0], [0, 0, 1], rb.rS + H, rb.rel);
GLX.setAll(pr, { u_model: m, u_center: rb.rel, u_R: rb.rS, u_H: H, u_atm: [...A.color, 1.0], u_k: A.k });
GLX.draw(this.meshes.mid, pr);
}
if (SOLAR) {
if (Settings.state.graphics.sunGlow) {
const s = this.sun, w = Math.max(V.len(sunR), 1e-9), px = this.pxPerUnit(w);
const minR = 26 / px; const size = Math.max(s.rS * 5.5, minR);
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'add' });
const pr = useP(P.corona);
GLX.setAll(pr, { u_center: sunR, u_size: size, u_scaleR: size / s.rS, u_int: lerp(1.4, 0.35, smoothstep(20, 260, s.proj.rpx)) * Gfx.flash() });
GLX.draw(this.meshes.quad, pr);
const pg = useP(P.glow), gsz = Math.max(s.rS * 14, 90 / px);
GLX.setAll(pg, { u_center: sunR, u_size: gsz, u_near: 0, u_color: [1.0, 0.62, 0.3], u_int: 0.06 * Math.min(1, (s.rS * px) / 12 + 0.25) * (1 - 0.9 * smoothstep(15, 120, s.proj.rpx)) * Gfx.flash(), u_sharp: 14 });
GLX.draw(this.meshes.quad, pg);
}
}
if (SOLAR) {
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const A = SCALES[ScaleState.from].glsl, B = SCALES[ScaleState.to].glsl;
const beltU = { u_days: this.jd - J2000, u_mapA: A, u_mapB: B, u_mapT: ScaleState.k(), u_off: sunR, u_dpr: GLX.dpr, u_ptK: ScaleState.dist(1) * 3 };
if (L.asteroids) { const pr = useP(P.belt); GLX.setAll(pr, beltU); GLX.set(pr, 'u_alpha', 0.55); GLX.draw(this.belts.ast, pr, Math.round(this.belts.ast.count * Gfx.particles())); }
if (L.kuiper) { const pr = useP(P.belt); GLX.setAll(pr, beltU); GLX.set(pr, 'u_ptK', ScaleState.dist(1) * 5); GLX.set(pr, 'u_alpha', 0.28); GLX.draw(this.belts.kui, pr, Math.round(this.belts.kui.count * Gfx.particles())); }
}
if (SOLAR && L.hz && vis1 > 0.003) {
if (!this.meshes.hz) {
const N = 512, pos = new Float32Array((N + 1) * 4), idx = new Uint16Array(N * 6);
for (let i = 0; i <= N; i++) { const a = i / N * Math.PI * 2; pos.set([0, a, 1, a], i * 4); }
for (let i = 0; i < N; i++) idx.set([i * 2, i * 2 + 1, i * 2 + 3, i * 2, i * 2 + 3, i * 2 + 2], i * 6);
this.meshes.hz = GLX.mesh({ a_pos: { data: pos, size: 2 } }, idx);
}
const H = HZ_MODELS[S.hzModel] || HZ_MODELS.conservador;
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const pr = useP(P.hz);
GLX.setAll(pr, { u_model: M4.translate(V.sub([0, 0, 0], cam)), u_rIn: ScaleState.dist(H.in), u_rOut: ScaleState.dist(H.out),
u_color: lin('#2f8f5f'), u_edge: lin('#8fd8aa'), u_alpha: 0.85 * vis1 });
GLX.draw(this.meshes.hz, pr);
}
if (SOLAR && typeof Planner !== 'undefined' && Planner.routeOn && Planner.route && vis1 > 0.003) {
const pts = Planner.scenePts(), n = pts ? Math.min(pts.length, 200) : 0;
if (n > 1) {
if (!this.meshes.route) { this._rtP = new Float32Array(600); this._rtF = new Float32Array(200); this.meshes.route = GLX.mesh({ a_pos: { data: this._rtP, size: 3 }, a_frac: { data: this._rtF, size: 1 } }, null, gl.LINE_STRIP, true); }
for (let i = 0; i < n; i++) { this._rtP[i * 3] = pts[i][0]; this._rtP[i * 3 + 1] = pts[i][1]; this._rtP[i * 3 + 2] = pts[i][2]; this._rtF[i] = i / (n - 1); }
GLX.update(this.meshes.route, 'a_pos', this._rtP, n); GLX.update(this.meshes.route, 'a_frac', this._rtF, n);
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const pr = useP(P.orbit);
GLX.setAll(pr, { u_model: M4.translate(V.sub([0, 0, 0], cam)), u_color: lin('#d39bff'), u_alpha: 0.95 * vis1, u_cur: (t * 0.12) % 1, u_fade: 0.7 });
GLX.draw(this.meshes.route, pr);
}
}
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const pO = useP(P.orbit);
const orbMode = S.mode === 'orbits';
if (vis1 > 0.003) for (const rb of this.rb) {
if (!rb.def.orbit) continue;                  // cuerpos centrales (Sol, Gargantúa): sin órbita
if (rb.isSun) continue;
const sel = UI.sel === rb || UI.tourRB === rb, hov = UI.hov === rb, type = rb.def.type;
let show = false, color = COLORS.orbitPlanet, alpha = 0.3;
if (type === 'craft') {
const o = rb.def.orbit;
if (!L.craft || rb.hidden || o.t === 'drift') continue;
if (o.t === 'leo') { const px = this.pxPerUnit(Math.max(V.dist(rb.parent.posS, cam), 1e-9)) * rb.moonDistMean; show = (L.moonOrbits && px > 14) || sel; alpha = 0.3 * smoothstep(14, 40, px); }
else show = L.orbits || sel;
color = COLORS.orbitCraft; alpha = o.t === 'leo' ? alpha : 0.32;
} else if (type === 'moon') {
if (!L.moons) continue;
const px = this.pxPerUnit(Math.max(V.dist(rb.parent.posS, cam), 1e-9)) * rb.moonDistMean;
show = (L.moonOrbits && px > 14) || sel; color = COLORS.orbitMoon; alpha = 0.3 * smoothstep(14, 40, px);
} else if (type === 'comet') { show = L.cometOrbits || orbMode || sel; color = COLORS.orbitComet; alpha = 0.32; }
else if (type === 'planet') { show = L.orbits || sel; alpha = 0.34; if (rb.def.fiction) { color = COLORS.orbitFiction; alpha = 0.26 * clamp((V.len(cam) / SYSTEMS.gargantua.H - 6) / 6, 0, 1); } }   // muy cerca de Gargantúa las órbitas se atenúan   // órbitas cinematográficas: otro color
else if (type === 'asteroid') { show = (L.orbits && orbMode) || sel; color = COLORS.orbitAst; alpha = 0.18; }
else { show = L.orbits || sel; color = COLORS.orbitDwarf; alpha = type === 'tno' ? 0.1 : 0.18; }
if (!show) continue;
if (orbMode) alpha = Math.min(1, alpha * 1.9);
if (hov) alpha = Math.max(alpha, 0.6);
if (sel) { color = COLORS.select; alpha = 0.9; }
let m, mesh;
if (type === 'craft' && rb.def.orbit.t === 'lpoint') { const Hh = rb.halo; m = M4.fromBasis(Hh.X, Hh.N, Hh.Z, 1, V.sub(V.add(rb.parent.posS, Hh.c), cam)); mesh = this.meshes.circle; }
else if (type === 'moon' || (type === 'craft' && rb.def.orbit.t === 'leo')) { const F = rb.oframe; m = M4.fromBasis(F.u, F.n, V.cross(F.u, F.n), rb.moonDistMean, V.sub(rb.parent.posS, cam)); mesh = this.meshes.circle; }
else { m = M4.translate(V.sub([0, 0, 0], cam)); mesh = rb.orbitMesh; }
GLX.setAll(pO, { u_model: m, u_color: lin(color), u_alpha: alpha * vis1, u_cur: rb.orbitFrac || 0, u_fade: orbMode ? 0.35 : 1 });
GLX.draw(mesh, pO);
}
if (SOLAR) {
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'add' });
for (const c of this.comets) this.drawComet(c, glob, useP);
if (Flight.on) Flight.drawFx(useP);
}
if (L.markers) {
const mb = this.markBuf; let n = 0;
for (const rb of this.rb) {
const show = rb.proj.on && UI.labelAllowed(rb, true);
const a = show ? (1 - smoothstep(1.5, 4, rb.proj.rpx)) * (rb.isSun ? 0 : 0.9) * vis1 : 0;
const rel = V.sub(rb.posS, cam);
mb.pos[n * 3] = rel[0]; mb.pos[n * 3 + 1] = rel[1]; mb.pos[n * 3 + 2] = rel[2];
const c = UI.sel === rb ? lin(COLORS.select) : rb.color;
mb.col[n * 4] = c[0] * 1.6; mb.col[n * 4 + 1] = c[1] * 1.6; mb.col[n * 4 + 2] = c[2] * 1.6; mb.col[n * 4 + 3] = a;
n++;
}
GLX.update(this.markers, 'a_pos', mb.pos); GLX.update(this.markers, 'a_col', mb.col);
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const pr = useP(P.mark); GLX.setAll(pr, { u_dpr: GLX.dpr, u_size: 5 });
GLX.draw(this.markers, pr);
}
GLX.endScene({ progs: this.post.progs, bloom: Settings.state.graphics.bloom ? lerp(0.7, 0.18, smoothstep(30, 260, this.sun.proj.on ? this.sun.proj.rpx : 0)) * Gfx.flash() : 0, threshold: GLX.hdr ? 1.1 : 0.27, exposure: 1.05, vignette: 0.35 });
},
drawComet(c, glob, useP) {
const rb = c.rb, cam = this.cam;
if (c.act < 0.01 || !rb.proj.on) return;
const r = c.r, h = rb.helio;
const anti = V.norm(h);
const h2 = Astro.helioPos(rb.def.orbit, this.jd + 0.5), vel = V.norm(V.sub(h2, h));
const LenAU = 0.12 + 0.55 * c.act;
const startS = rb.posS;
const endS = ScaleState.mapVec(Astro.eclToScene(V.add(h, V.scale(anti, LenAU))));
const lenS = V.dist(startS, endS);
const N = c.N;
const fill = (mesh, kind) => {
const pos = new Float32Array((N + 1) * 6), uv = new Float32Array((N + 1) * 4);
let prev = null;
for (let i = 0; i <= N; i++) {
const s = i / N;
let dir = anti;
if (kind === 'dust') dir = V.norm(V.sub(anti, V.scale(vel, 0.55 * s)));
const pAU = V.add(h, V.scale(dir, LenAU * s * (kind === 'dust' ? 0.8 : 1)));
const pS = ScaleState.mapVec(Astro.eclToScene(pAU));
const rel = V.sub(pS, cam);
const tan = prev ? V.norm(V.sub(rel, prev)) : V.norm(V.sub(ScaleState.mapVec(Astro.eclToScene(V.add(h, V.scale(dir, LenAU * 0.02)))), startS));
prev = rel;
let side = V.cross(tan, V.norm(rel));
const sl = V.len(side); side = sl > 1e-9 ? V.scale(side, 1 / sl) : [0, 1, 0];
const wdt = lenS * (kind === 'dust' ? 0.02 + 0.2 * s : 0.008 + 0.06 * s);
for (let k = 0; k < 2; k++) {
const sg = k === 0 ? -1 : 1, j = (i * 2 + k);
pos[j * 3] = rel[0] + side[0] * wdt * sg; pos[j * 3 + 1] = rel[1] + side[1] * wdt * sg; pos[j * 3 + 2] = rel[2] + side[2] * wdt * sg;
uv[j * 2] = s; uv[j * 2 + 1] = sg;
}
}
GLX.update(mesh, 'a_pos', pos); GLX.update(mesh, 'a_uv', uv);
};
fill(c.ion, 'ion'); fill(c.dust, 'dust');
const w = Math.max(V.dist(startS, cam), 1e-9), px = this.pxPerUnit(w);
const pr = useP(this.P.tail);
const nearT = lenS * 0.025;
GLX.setAll(pr, { u_color: [0.35, 0.6, 1.0], u_int: 0.4 * c.act, u_seed: c.seed, u_near: nearT });
GLX.draw(c.ion, pr);
GLX.setAll(pr, { u_color: [1.0, 0.86, 0.66], u_int: 0.3 * c.act, u_seed: c.seed + 1.0, u_near: nearT });
GLX.draw(c.dust, pr);
const comaS = Math.min(lenS * 0.04, rb.rS * (25 + 120 * Math.min(c.act, 1)));
const size = Math.max(comaS, 6 / px);
const sizePx = size * px;
const k = Math.min(1, c.act) * (1 - smoothstep(40, 220, sizePx));
if (k < 0.003) return;
const pg = useP(this.P.glow);
GLX.setAll(pg, { u_center: V.sub(startS, cam), u_size: size, u_color: [0.55, 1.0, 0.8], u_int: 0.34 * k, u_sharp: 12, u_near: size });
GLX.draw(this.meshes.quad, pg);
},
};
const EarthMap = {
build() {
const W = 2048, H = 1024;
const mk = (w, h) => { const c = document.createElement('canvas'); c.width = w; c.height = h; return c; };
const toXY = (lon, lat, w, h) => [(lon + 180) / 360 * w, (90 - lat) / 180 * h];
const poly = (ctx, str, w, h) => {
ctx.beginPath();
str.trim().split(/\s+/).forEach((pt, i) => { const [lo, la] = pt.split(',').map(Number); const [x, y] = toXY(lo, la, w, h); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
ctx.closePath(); ctx.fill();
};
const land = mk(W, H), lc = land.getContext('2d');
lc.fillStyle = '#000'; lc.fillRect(0, 0, W, H);
lc.fillStyle = '#fff'; EARTH_LAND.forEach(s => poly(lc, s, W, H));
lc.fillStyle = '#000'; EARTH_WATER.forEach(s => poly(lc, s, W, H));
const blurVia = (src, sw, sh) => { const sm = mk(sw, sh), sc = sm.getContext('2d'); sc.imageSmoothingEnabled = true; sc.drawImage(src, 0, 0, sw, sh); const big = mk(W, H), bc = big.getContext('2d'); bc.imageSmoothingEnabled = true; bc.drawImage(sm, 0, 0, W, H); return bc.getImageData(0, 0, W, H).data; };
const inland = blurVia(land, 96, 48);
const mtn = mk(W, H), mc = mtn.getContext('2d');
mc.fillStyle = '#000'; mc.fillRect(0, 0, W, H);
mc.strokeStyle = '#fff'; mc.lineWidth = 11; mc.lineCap = 'round'; mc.lineJoin = 'round';
EARTH_RANGES.forEach(s => { mc.beginPath(); s.split(/\s+/).forEach((pt, i) => { const [lo, la] = pt.split(',').map(Number); const [x, y] = toXY(lo, la, W, H); if (i) mc.lineTo(x, y); else mc.moveTo(x, y); }); mc.stroke(); });
const mBlur = blurVia(mtn, 256, 128);
const ld = lc.getImageData(0, 0, W, H);
const out = mk(W, H), oc = out.getContext('2d'), od = oc.createImageData(W, H);
for (let i = 0; i < W * H; i++) {
const L = ld.data[i * 4];
od.data[i * 4] = L;
od.data[i * 4 + 1] = Math.min(255, inland[i * 4] * 1.15);
od.data[i * 4 + 2] = L > 128 ? Math.min(255, mBlur[i * 4] * 1.6) : 0;
od.data[i * 4 + 3] = 255;
}
oc.putImageData(od, 0, 0);
const fin = mk(W, H), fc = fin.getContext('2d');
fc.drawImage(blurSmall(out, W, H), 0, 0);
return fin;
function blurSmall(c, w, h) { const s = mk(w / 2, h / 2); s.getContext('2d').drawImage(c, 0, 0, w / 2, h / 2); const b = mk(w, h), bx = b.getContext('2d'); bx.imageSmoothingEnabled = true; bx.drawImage(s, 0, 0, w, h); return b; }
},
};
const CraftModels = (() => {
const C = {
gold: [0.86, 0.62, 0.22, 0.9], silver: [0.5, 0.52, 0.56, 0.5], white: [0.72, 0.73, 0.75, 0.25], dark: [0.12, 0.12, 0.13, 0.2],
panel: [0.08, 0.12, 0.3, 0.7], panelGold: [0.6, 0.38, 0.14, 0.6], kapton: [0.62, 0.55, 0.68, 0.6], mirror: [0.95, 0.72, 0.25, 1.2], black: [0.03, 0.03, 0.035, 0.1],
};
function Builder() { this.p = []; this.n = []; this.c = []; }
Builder.prototype.tri = function (a, b, c, col) {
let n = V.norm(V.cross(V.sub(b, a), V.sub(c, a)));
for (const v of [a, b, c]) { this.p.push(v[0], v[1], v[2]); this.n.push(n[0], n[1], n[2]); this.c.push(col[0], col[1], col[2], col[3]); }
};
Builder.prototype.quad = function (a, b, c, d, col) { this.tri(a, b, c, col); this.tri(a, c, d, col); };
Builder.prototype.box = function (ctr, size, col, ax) {
ax = ax || [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const h = size.map(s => s / 2), P = (sx, sy, sz) => V.add(ctr, V.add(V.add(V.scale(ax[0], sx * h[0]), V.scale(ax[1], sy * h[1])), V.scale(ax[2], sz * h[2])));
const v = [P(-1, -1, -1), P(1, -1, -1), P(1, 1, -1), P(-1, 1, -1), P(-1, -1, 1), P(1, -1, 1), P(1, 1, 1), P(-1, 1, 1)];
const f = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]];
f.forEach(([a, b, c, d]) => this.quad(v[a], v[b], v[c], v[d], col));
};
Builder.prototype.cyl = function (ctr, axis, r, len, seg, col, capCol) {
axis = V.norm(axis);
let u = V.norm(V.cross(axis, Math.abs(axis[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0])), w = V.cross(axis, u);
const a0 = V.add(ctr, V.scale(axis, -len / 2)), a1 = V.add(ctr, V.scale(axis, len / 2));
const ring = (o, k) => V.add(o, V.add(V.scale(u, Math.cos(TAU * k / seg) * r), V.scale(w, Math.sin(TAU * k / seg) * r)));
for (let k = 0; k < seg; k++) {
const p0 = ring(a0, k), p1 = ring(a0, k + 1), q0 = ring(a1, k), q1 = ring(a1, k + 1);
this.quad(p0, p1, q1, q0, col);
this.tri(a1, q0, q1, capCol || col); this.tri(a0, p1, p0, capCol || col);
}
};
Builder.prototype.dish = function (ctr, axis, r, depth, seg, col, backCol) {
axis = V.norm(axis);
let u = V.norm(V.cross(axis, Math.abs(axis[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0])), w = V.cross(axis, u);
const apex = V.add(ctr, V.scale(axis, -depth));
const rim = k => V.add(ctr, V.add(V.scale(u, Math.cos(TAU * k / seg) * r), V.scale(w, Math.sin(TAU * k / seg) * r)));
for (let k = 0; k < seg; k++) { this.tri(apex, rim(k + 1), rim(k), col); this.tri(apex, rim(k), rim(k + 1), backCol || col); }
};
Builder.prototype.poly = function (ctr, u, w, pts, col, backCol, thick) {
const n = V.norm(V.cross(u, w)), t = (thick || 0.02) / 2;
const P = (q, s) => V.add(V.add(ctr, V.scale(n, s * t)), V.add(V.scale(u, q[0]), V.scale(w, q[1])));
for (let k = 1; k < pts.length - 1; k++) {
this.tri(P(pts[0], 1), P(pts[k], 1), P(pts[k + 1], 1), col);
this.tri(P(pts[0], -1), P(pts[k + 1], -1), P(pts[k], -1), backCol || col);
}
};
const hex = (r) => [0, 1, 2, 3, 4, 5].map(k => [Math.cos(k * Math.PI / 3 + Math.PI / 6) * r, Math.sin(k * Math.PI / 3 + Math.PI / 6) * r]);
const MODELS = {
jwst(B) {
const shield = [[10.6, 0], [5.2, 7.1], [-5.2, 7.1], [-10.6, 0], [-5.2, -7.1], [5.2, -7.1]];
for (let i = 0; i < 5; i++) {
const s = 1 - i * 0.035;
B.poly([0, i * 0.35, 0], [1, 0, 0], [0, 0, -1], shield.map(q => [q[0] * s, q[1] * s]), i === 0 ? C.silver : C.kapton, C.silver, 0.04);
}
const size = 0.762, u = V.norm([1, 0, 0]), w = V.norm([0, Math.cos(0.35), Math.sin(0.35)]), ctr = [0, 6.0, 1.2];
for (let q = -2; q <= 2; q++) for (let r = -2; r <= 2; r++) {
const s = -q - r, d = Math.max(Math.abs(q), Math.abs(r), Math.abs(s));
if (d === 0 || d > 2) continue;
if (d === 2 && (Math.abs(q) === 2 && Math.abs(r) === 0 || 0) && false) continue;
const x = size * Math.sqrt(3) * (q + r / 2), y = size * 1.5 * r;
const n = V.norm(V.cross(u, w));
B.poly(V.add(ctr, V.add(V.scale(u, x), V.scale(w, y))), u, w, hex(size * 0.96), C.mirror, C.black, 0.08);
B.p.length; void n;
}
const nrm = V.norm(V.cross(u, w)), sec = V.add(ctr, V.scale(nrm, -7.0));
B.cyl(sec, nrm, 0.37, 0.15, 12, C.mirror, C.black);
[[-2.2, -1.6], [2.2, -1.6], [0, 2.6]].forEach(([x, y]) => {
const base = V.add(ctr, V.add(V.scale(u, x), V.scale(w, y)));
const mid = V.scale(V.add(base, sec), 0.5), dir = V.sub(sec, base);
B.cyl(mid, dir, 0.06, V.len(dir), 5, C.dark);
});
B.cyl([0, 3.0, 0.6], [0, 1, 0], 0.35, 5.6, 6, C.dark);
B.box([0, -1.4, 0], [3.0, 1.6, 3.0], C.dark);
B.box([0, -1.4, -3.4], [1.2, 0.08, 3.8], C.panel);
},
hubble(B) {
B.cyl([0, 0, 0], [1, 0, 0], 2.1, 13.2, 20, C.silver, C.dark);
B.cyl([-7.4, 0, 0], [1, 0, 0], 2.1, 1.6, 20, C.dark);
[-1, 1].forEach(sg => {
B.box([0.5, 0, sg * 5.6], [2.6, 0.06, 7.2], C.panel);
B.cyl([0.5, 0, sg * 2.4], [0, 0, 1], 0.08, 2.0, 5, C.silver);
});
B.cyl([3.0, 2.6, 0], [0, 1, 0], 0.06, 1.2, 5, C.white); B.dish([3.0, 3.3, 0], [0, 1, 0], 0.65, 0.2, 14, C.white);
},
iss(B) {
B.box([0, 0, 0], [109, 1.8, 1.8], C.silver);
[-1, 1].forEach(sx => [38, 50].forEach(x => [-1, 1].forEach(sz => {
B.box([sx * x, 0, sz * 19], [11.6, 0.08, 34], C.panelGold);
B.box([sx * x, 0, sz * 19], [0.5, 0.5, 34], C.dark);
})));
[-1, 1].forEach(sx => B.box([sx * 20, -6, 0], [3.2, 0.1, 12], C.white));
B.cyl([0, -4.5, -4], [0, 0, 1], 2.1, 46, 14, C.white);
B.cyl([0, -4.5, 13], [1, 0, 0], 2.1, 18, 14, C.white);
B.cyl([0, -4.5, -24], [1, 0, 0], 2.1, 14, 14, C.white);
[-1, 1].forEach(s => B.box([0, -4.5, -38 + s * 0.0], [s * 0 + 3, 0.06, 20], C.panel));
},
tiangong(B) {
B.cyl([0, 0, 0], [0, 0, 1], 2.1, 16.6, 14, C.white);
[-1, 1].forEach(s => {
B.cyl([s * 9.0, 0, 8.6], [1, 0, 0], 2.1, 17.9, 14, C.white);
B.box([s * 13, 0, 8.6], [2.5, 0.06, 27], C.panel);
});
B.box([0, 0, -6], [27, 0.06, 2.3], C.panel);
},
soho(B) {
B.box([0, 0, 0], [4.3, 3.7, 2.7], C.gold);
B.box([0, 2.0, 0], [3.0, 0.4, 2.0], C.white);
[-1, 1].forEach(s => B.box([s * 4.45, 0, 0], [4.6, 0.06, 1.3], C.panel));
B.dish([0, -1.85, 1.2], [0, -1, 0.5], 0.6, 0.18, 12, C.white);
},
euclid(B) {
B.cyl([0, 0, 0], [0, 0, 1], 1.4, 3.6, 8, C.silver, C.dark);
B.box([0, 1.75, 0], [3.4, 0.08, 4.7], C.panel);
B.box([0, -0.9, -1.9], [2.6, 1.2, 0.9], C.gold);
},
voyager(B) {
B.dish([0, 1.0, 0], [0, 1, 0], 1.85, 0.45, 24, C.white, C.silver);
B.cyl([0, 0.15, 0], [0, 1, 0], 0.9, 0.5, 10, C.dark);
B.cyl([0, 0.55, 0], [0, 1, 0], 0.1, 0.5, 6, C.dark);
B.cyl([-6.5, 0, 0], [1, 0, 0], 0.03, 13, 4, C.silver);
B.cyl([1.8, 0.1, 1.0], [1, 0, 0.6], 0.04, 2.5, 4, C.silver);
[0, 1, 2].forEach(k => B.cyl([2.5 + k * 0.42, 0.1, 1.45 + k * 0.25], [1, 0, 0.6], 0.2, 0.4, 8, C.black));
B.cyl([0.2, 0.1, -2.2], [0, 0, -1], 0.05, 2.4, 4, C.silver);
B.box([0.2, 0.1, -3.4], [0.4, 0.4, 0.5], C.dark);
},
newhorizons(B) {
B.poly([0, 0, 0], [1, 0, 0], [0, 0, -1], [[1.3, -0.75], [0, 1.5], [-1.3, -0.75]], C.gold, C.gold, 0.7);
B.dish([0, 0.65, 0], [0, 1, 0], 1.05, 0.25, 20, C.white, C.silver);
B.cyl([1.6, 0, -0.9], [1, 0, -0.5], 0.22, 1.2, 10, C.black);
},
parker(B) {
B.cyl([0, 1.2, 0], [0, 1, 0], 1.15, 0.12, 24, C.dark, C.white);
B.cyl([0, 0.15, 0], [0, 1, 0], 0.5, 1.6, 6, C.gold);
[-1, 1].forEach(s => B.box([s * 1.0, -0.2, 0], [0.9, 0.04, 0.45], C.panel));
B.cyl([0, -1.1, 0], [0, 1, 0], 0.04, 1.4, 4, C.silver);
},
};
const cache = {};
return {
mesh(name) {
if (cache[name]) return cache[name];
if (!MODELS[name]) return null;          // sin versión procedural (naves con modelo importado únicamente)
const B = new Builder(); MODELS[name](B, C);
let R = 0; for (let i = 0; i < B.p.length; i += 3) R = Math.max(R, Math.hypot(B.p[i], B.p[i + 1], B.p[i + 2]));
const pos = new Float32Array(B.p.map(v => v / R));
const m = GLX.mesh({ a_pos: { data: pos, size: 3 }, a_nrm: { data: new Float32Array(B.n), size: 3 }, a_col: { data: new Float32Array(B.c), size: 4 } });
m.R = R;
return (cache[name] = m);
},
define(name, fn) { MODELS[name] = fn; },
setMesh(name, mesh) { if (MODELS[name] && !MODELS[name + '-lite']) MODELS[name + '-lite'] = MODELS[name]; cache[name] = mesh; },
hex,
};
})();
COLORS.orbitFiction = '#c9a46e';
Object.assign(World, {
system: 'solar', bhReveal: 1,
buildSystems() {
this.systems = { solar: { rb: this.rb, sun: this.sun } };
const G = SYSTEMS.gargantua, list = [];
GARG_BODIES.forEach(def => { const rb = this.makeBody(def); rb.R = def.visR; rb.extent = def.extent || null; list.push(rb); this.byId[rb.id] = rb; });
list.forEach(rb => { if (rb.def.parent) { rb.parent = this.byId[rb.def.parent]; rb.parent.children.push(rb); } });
const bh = this.byId.gargantua; bh.isBH = true; bh.rS = G.H;
list.forEach(rb => {
const o = rb.def.orbit; if (!o) return;
const N = 360, pos = new Float32Array((N + 1) * 3), fr = new Float32Array(N + 1);
for (let i = 0; i <= N; i++) { const a = i / N * Math.PI * 2, p = this.cinePos(o, a); pos.set(p, i * 3); fr[i] = i / N; }
rb.orbitMesh = GLX.mesh({ a_pos: { data: pos, size: 3 }, a_frac: { data: fr, size: 1 } }, null, GLX.gl.LINE_STRIP);
});
this.systems.gargantua = { rb: list, sun: bh };
},
cinePos(o, a) { const x = o.a * Math.cos(a), z = -o.a * Math.sin(a); return [x, z * Math.sin(o.i || 0), z * Math.cos(o.i || 0)]; },
setSystem(id) {
const s = this.systems[id]; if (!s || id === this.system) return;
this.system = id; this.rb = s.rb; this.sun = s.sun;
this.rb.forEach(rb => { rb.hover = 0; rb.drawOn = false; });
for (const k in this.systems) this.systems[k].rb.forEach(rb => { const L = UI.labels[rb.id]; if (L) { L.el.style.display = k === id ? '' : 'none'; L.shown = false; } });
},
updateSystem(jd, dtReal, speed) {
this.jd = jd; const d = jd - J2000;
for (const rb of this.rb) {
const o = rb.def.orbit; rb.rS = rb.def.visR;
if (o && o.t === 'cine') { const a = (o.phase + d / o.P) * Math.PI * 2; rb.posS = this.cinePos(o, a); rb.orbitFrac = ((a / (Math.PI * 2)) % 1 + 1) % 1; }
else rb.posS = [0, 0, 0];
rb.helio = V.scale(rb.posS, 1 / 1000); rb.km = V.scale(rb.posS, 1000);
rb.rotPhase = (rb.def.W0 || 0) + rb.rate * d; rb.basis = Astro.basisFromFrame(rb.frame, rb.rotPhase);
}
const bh = this.sun, cam = this.cam || [0, 0, 0], rel = V.sub(bh.posS, cam), dB = V.len(rel), sh = Math.atan(bh.rS * 2.6 / Math.max(dB, 1e-6));
for (const rb of this.rb) { if (rb === bh) continue; const rp = V.sub(rb.posS, cam), dp = V.len(rp); rb.occBH = dp > dB && Math.acos(clamp(V.dot(rel, rp) / (dB * dp), -1, 1)) < sh; }
},
drawBH(t, useP) {
const G = SYSTEMS.gargantua, bh = this.sun, pr = useP(this.P.bh), v = this.view, f = 1 / this.proj[5];
const preset = Settings.state.graphics.preset; let steps = { low: 60, medium: 100, high: 150, ultra: 220 }[preset] || 130;
const dist = Math.max(V.len(V.sub(this.cam, bh.posS)), 1e-6), cover = clamp((G.disk.out * 1.7 * G.H / dist) / f, 0, 3);
steps = Math.round(steps * clamp(1.35 - cover * 0.45, 0.45, 1));
GLX.state({ depthTest: false, depthWrite: false, cull: 'none', blend: 'premul' });
GLX.setAll(pr, { u_camR: [v[0], v[4], v[8]], u_camU: [v[1], v[5], v[9]], u_camF: [-v[2], -v[6], -v[10]], u_tan: [1 / this.proj[0], 1 / this.proj[5]], u_off: [this.proj[8], this.proj[9]],
u_cam: V.scale(V.sub(this.cam, bh.posS), 1 / G.H), u_steps: steps, u_rin: G.disk.in, u_rout: G.disk.out, u_reveal: this.bhReveal,
u_heat: 1 + (this.bhBoost || 0), u_spin: UI.reducedMotion ? 0.12 : 0.45 });
GLX.gl.bindVertexArray(GLX.fsTri); GLX.gl.drawArrays(GLX.gl.TRIANGLES, 0, 3);
},
});
const Cam = {
target: [0, 0, 0], dist: 560, dDist: 560, az: 0.9, el: 0.36, fov: 42 * DEG,
vAz: 0, vEl: 0, focus: null, follow: true, fly: null, look: null, auto: 0,
pos: [0, 0, 0],
offsetDir() { const ce = Math.cos(this.el); return [ce * Math.sin(this.az), Math.sin(this.el), ce * Math.cos(this.az)]; },
minDist() { const f = this.focus; if (!f) return 1e-6; if (f.extent) return f.rS * f.extent * 1.18;   // forma importada: según su radio máximo
return f.rS * (f.irregular ? 1.6 : f.isCraft ? 1.3 : 1.12); },
maxDist() { return World.system !== 'solar' ? SYSTEMS[World.system].maxDist : ScaleState.get('maxDist'); },   // cada sistema tiene su propio límite
frameDist(rb, close) {
let k = (close ? 2.4 : 4.2) * (rb.extent || 1);
if (rb.def.vis.rings && rb.def.vis.rings.type === 'saturn') k *= 1.7;
if (rb.def.model === 'voyager' && World.hd.voyager) k *= 0.55;
if (rb.isSun) k = close ? 3 : 5;
if (innerWidth < 760) k *= 1.5 * Math.max(1, (innerHeight / innerWidth) * 0.75);
return rb.rS * k;
},
arrivalAngles(rb) {
if (rb.isSun) return { az: this.az, el: 0.3 };
if (rb.def.model === 'jwst') { const s = V.norm(V.scale(rb.posS, -1)); return { az: Math.atan2(s[0], s[2]) + 2.1, el: 0.55 }; }
const s = V.norm(V.scale(rb.posS, -1));
return { az: Math.atan2(s[0], s[2]) + 0.75, el: 0.22 + Math.asin(clamp(s[1], -1, 1)) * 0.5 };
},
travel(rb, opts) {
opts = Object.assign({}, opts || {});
if (World.system === 'gargantua') opts.cine = true;      // Sistema Gargantúa: transiciones más lentas y dramáticas
const d1 = opts.dist || this.frameDist(rb, opts.close);
const ang = opts.keepAngles ? { az: this.az, el: this.el } : this.arrivalAngles(rb);
const travel = V.dist(this.target, rb.posS);
const d0 = this.dist;
const ratio = Math.log(1 + travel / (d0 + d1));
const reduce = UI.reducedMotion;
let dur = clamp(1.3 + 0.55 * ratio + 0.12 * Math.abs(Math.log(d1 / d0)), 1.2, 4.2);
if (opts.dur) dur = opts.dur;
if (opts.cine) dur = clamp(dur * 1.8, 3.4, 7.5);
if (reduce) dur = 0.45;
let daz = wrapPi(ang.az - this.az);
const newFocus = this.focus !== rb;
this.fly = { t: 0, dur, from: { target: this.target.slice(), dist: d0, az: this.az, el: this.el }, body: rb, d1, az1: this.az + daz, el1: ang.el, arc: reduce ? 0 : clamp(ratio, 0, 3.5) * (opts.cine ? 0.7 : 0.55), cine: !!opts.cine };
this.focus = rb; this.follow = true; this.vAz = this.vEl = 0; this.look = null;
if (newFocus || ratio > 0.6) { this.fly.sfx = true; SFX.navStart(this.fly, ratio / 3.5); } else SFX.navCancel();
},
travelTo(target, dist, el) {
const reduce = UI.reducedMotion;
this.fly = { t: 0, dur: reduce ? 0.45 : 2.2, from: { target: this.target.slice(), dist: this.dist, az: this.az, el: this.el }, point: target.slice(), d1: dist, az1: this.az + 0.25, el1: el, arc: 0, sfx: true };
this.focus = null; this.follow = false; this.look = null;
SFX.navStart(this.fly, 0.45);
},
lookDir(D) { // orientar la cámara para mirar en la dirección D
const az = Math.atan2(-D[0], -D[2]), el = Math.asin(clamp(-D[1], -1, 1));
this.look = { t: 0, dur: UI.reducedMotion ? 0.3 : 1.6, az0: this.az, el0: this.el, az1: this.az + wrapPi(az - this.az), el1: clamp(el, -1.5, 1.5) };
},
center() { if (this.focus) this.travel(this.focus, { dist: clamp(this.dist, this.minDist() * 1.05, this.maxDist()), keepAngles: true, dur: 1.0 }); },
zoom(f) { this.dDist = clamp(this.dDist * f, this.minDist(), this.maxDist()); if (this.fly && this.fly.sfx) SFX.navCancel(); this.fly = null; },
cancel() { if (this.fly && this.fly.sfx) SFX.navCancel(); this.fly = null; this.look = null; },
update(dt) {
if (this.fly) {
const F = this.fly; F.t += dt;
const u = clamp(F.t / F.dur, 0, 1);
const k = F.cine ? (u < 0.5 ? 16 * u ** 5 : 1 - Math.pow(-2 * u + 2, 5) / 2) : easeInOut(u);   // Modo Cine: arranque y llegada muy suaves
const dest = F.body ? F.body.posS : F.point;
this.target = V.lerp(F.from.target, dest, k);
this.dist = Math.exp(lerp(Math.log(F.from.dist), Math.log(F.d1), k) + F.arc * Math.sin(Math.PI * k));
this.dDist = this.dist;
this.az = lerp(F.from.az, F.az1, k); this.el = lerp(F.from.el, F.el1, k);
if (F.sfx) SFX.navUpdate(F.t / F.dur);
if (F.t >= F.dur) { this.fly = null; if (F.body) { this.focus = F.body; this.follow = true; } if (F.sfx) SFX.navArrive(); }
} else {
if (this.focus && this.follow) this.target = this.focus.posS.slice();
const md = this.minDist();
if (this.dDist < md) this.dDist = md;
this.dist = Math.exp(lerp(Math.log(Math.max(this.dist, 1e-9)), Math.log(this.dDist), 1 - Math.exp(-dt * 9)));
if (this.look) {
const L = this.look; L.t += dt; const k = easeInOut(clamp(L.t / L.dur, 0, 1));
this.az = lerp(L.az0, L.az1, k); this.el = lerp(L.el0, L.el1, k);
if (L.t >= L.dur) this.look = null;
}
if (!Input.dragging) {
this.az += this.vAz * dt; this.el += this.vEl * dt;
const damp = Math.exp(-dt * 4); this.vAz *= damp; this.vEl *= damp;
}
if (this.auto) this.az += this.auto * dt;
}
this.el = clamp(this.el, -1.53, 1.53);
const o = this.offsetDir();
this.pos = V.add(this.target, V.scale(o, this.dist));
World.cam = this.pos;
const f = V.scale(o, -1);
World.view = M4.viewRot(f, [0, 1, 0], World.view);
const aspect = GLX.W / GLX.H;
World.proj = M4.perspInf(this.fov, aspect, clamp(this.dist * 2e-4, 1e-9, 1), World.proj);
const W = innerWidth, H = innerHeight, info = UI.sel && document.getElementById('info').classList.contains('open') && S.mode !== 'cine' && !UI.tour;
let sx = 0, sy = 0;
if (W >= 760) {
const tourSide = UI.tour && !UI.tour.done && W >= 900;
const R = info ? 396 : tourSide ? 408 : 0, L = (UI.section && W > 980) ? 414 : 0;
sx = (R - L) / W;
} else if (info) sy = -(0.56 * H) / H;
if (UI.tour && !UI.tour.done && W < 900) { const tb = document.getElementById('tourbar'); if (tb) sy = -Math.min(tb.offsetHeight + 40, H * 0.6) / H; }
const k = 1 - Math.exp(-dt * 5);
this.sx = lerp(this.sx || 0, sx, k); this.sy = lerp(this.sy || 0, sy, k);
this.isy = lerp(this.isy || 0, this.introSY || 0, k);              // composición de presentaciones (p. ej., Gargantúa en el tercio superior)
World.proj[8] = this.sx; World.proj[9] = this.sy + this.isy;
},
};
const Input = {
pointers: new Map(), dragging: false, mode: null, last: null, down: null, pinch: null,
init(cv) {
cv.addEventListener('contextmenu', e => e.preventDefault());
cv.addEventListener('pointerdown', e => { if (!Flight.on) this.onDown(e); });
cv.addEventListener('pointermove', e => { if (!Flight.on) this.onMove(e); });
cv.addEventListener('pointerup', e => { if (!Flight.on) this.onUp(e); });
cv.addEventListener('pointercancel', e => { if (!Flight.on) this.onUp(e); });
cv.addEventListener('pointerleave', () => { if (!this.dragging && !Flight.on) UI.setHover(null); });
cv.addEventListener('wheel', e => {
e.preventDefault();
if (Flight.on) return;
let dy = e.deltaY; if (e.deltaMode === 1) dy *= 16; else if (e.deltaMode === 2) dy *= 400;
Cam.zoom(Math.exp(clamp(dy, -300, 300) * 0.0013));
UI.camTouched();
}, { passive: false });
cv.addEventListener('dblclick', e => {
if (Flight.on) return;
const rb = World.pick(e.offsetX, e.offsetY);
if (rb) { UI.select(rb, { fly: true, close: true }); }
});
},
onDown(e) {
e.target.setPointerCapture(e.pointerId);
this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
this.down = { x: e.clientX, y: e.clientY, t: performance.now(), moved: 0 };
this.last = { x: e.clientX, y: e.clientY };
this.mode = (e.button === 2 || e.shiftKey || e.ctrlKey || e.metaKey) ? 'pan' : 'rot';
if (this.pointers.size === 2) { this.mode = 'pinch'; this.pinch = this.pinchState(); }
this.dragging = true; Cam.vAz = Cam.vEl = 0;
UI.camTouched();
},
pinchState() { const p = [...this.pointers.values()]; return { d: Math.hypot(p[0].x - p[1].x, p[0].y - p[1].y), cx: (p[0].x + p[1].x) / 2, cy: (p[0].y + p[1].y) / 2 }; },
onMove(e) {
if (!this.dragging) { const r = e.target.getBoundingClientRect(); UI.setHover(World.pick(e.clientX - r.left, e.clientY - r.top)); return; }
if (!this.pointers.has(e.pointerId)) return;
this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
const dx = e.clientX - this.last.x, dy = e.clientY - this.last.y;
this.last = { x: e.clientX, y: e.clientY };
this.down.moved += Math.abs(dx) + Math.abs(dy);
if (this.mode === 'pinch' && this.pointers.size >= 2) {
const s = this.pinchState();
if (this.pinch.d > 0 && s.d > 0) Cam.zoom(this.pinch.d / s.d);
this.pan(s.cx - this.pinch.cx, s.cy - this.pinch.cy);
this.pinch = s; return;
}
if (this.down.moved < 3) return;
Cam.cancel();
if (this.mode === 'pan') this.pan(dx, dy);
else {
const k = 0.0055;
Cam.az -= dx * k; Cam.el += dy * k;
Cam.vAz = -dx * k / 0.016 * 0.35; Cam.vEl = dy * k / 0.016 * 0.35;
}
},
pan(dx, dy) {
if (Math.abs(dx) + Math.abs(dy) < 0.01) return;
const v = World.view, H = GLX.H / GLX.dpr;
const s = Cam.dist * 2 * Math.tan(Cam.fov / 2) / H;
const right = [v[0], v[4], v[8]], up = [v[1], v[5], v[9]];
Cam.target = V.add(Cam.target, V.add(V.scale(right, -dx * s), V.scale(up, dy * s)));
Cam.follow = false; Cam.fly = null;
UI.refreshCamButtons();
},
onUp(e) {
this.pointers.delete(e.pointerId);
if (this.pointers.size === 1 && this.mode === 'pinch') { this.mode = 'rot'; const p = [...this.pointers.values()][0]; this.last = { x: p.x, y: p.y }; return; }
if (this.pointers.size > 0) return;
const d = this.down;
this.dragging = false;
if (d && d.moved < 6 && performance.now() - d.t < 450 && e.type === 'pointerup') {
const r = e.target.getBoundingClientRect();
const rb = World.pick(e.clientX - r.left, e.clientY - r.top);
if (rb) UI.select(rb, { fly: true });
}
},
};
World.pick = function (x, y) {
let best = null, bestScore = Infinity;
for (const rb of this.rb) {
const p = rb.proj;
if (!p.on || !UI.labelAllowed(rb, true)) continue;
const d = Math.hypot(p.x - x, p.y - y), R = Math.max(p.rpx, 9);
if (d > R + 4) continue;
const inside = d < p.rpx;
const score = inside ? -1e6 - (1e3 / Math.max(p.w, 1e-9)) * 0 + p.w * 1e-3 : d + (rb.def.type === 'moon' ? 3 : 0);
if (score < bestScore) { bestScore = score; best = rb; }
}
return best;
};
const Time = {
jd: Astro.jdFromDate(new Date()),
speeds: [1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7],
idx: 0, paused: false, live: true,
get speed() { return this.paused ? 0 : this.speeds[this.idx]; },
step(dt) {
if (this.live && !this.paused && this.idx === 0) { this.jd = Astro.jdFromDate(new Date()); return; }
if (!this.paused) this.jd += dt * this.speeds[this.idx] / 86400;
},
goLive() { this.jd = Astro.jdFromDate(new Date()); this.idx = 0; this.paused = false; this.live = true; },
startState() {
const sim = Settings.state.simulation;
this.goLive();
if (sim.startDate === 'last' && sim.lastJD) { this.jd = sim.lastJD; this.live = false; }
if (sim.timeSpeed === 'paused') this.paused = true;
else if (sim.timeSpeed !== 'real') { this.idx = { x10: 1, x100: 2, x1000: 3 }[sim.timeSpeed] || 0; this.live = false; }
},
label() {
const s = this.speeds[this.idx];
return s === 1 ? 'Tiempo real' : '×' + s.toLocaleString(I18N.loc()).replace(/,/g, ' ');
},
rateText() {
const s = this.speeds[this.idx];
if (s < 60) return s + ' s por segundo';
if (s < 3600) return Math.round(s / 60) + ' min por segundo';
if (s < 86400) return (s / 3600).toFixed(1).replace('.0', '') + ' h por segundo';
return (s / 86400).toFixed(s < 864000 ? 1 : 0).replace('.0', '') + ' días por segundo';
},
};
const ICON = (() => {
const s = d => `<svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${d}</svg>`;
return {
explorar: s('<circle cx="12" cy="12" r="9"/><path d="M15.5 8.5l-2 5-5 2 2-5z"/>'),
sistema: s('<circle cx="12" cy="12" r="2.5"/><ellipse cx="12" cy="12" rx="9.5" ry="4.5"/><circle cx="20.5" cy="11" r="1"/>'),
planetas: s('<circle cx="12" cy="12" r="5"/><path d="M3.5 15.5c2-1 15-6 17-8.5"/>'),
lunas: s('<path d="M15.5 3.5a8.5 8.5 0 1 0 5 13 7 7 0 0 1-5-13z"/>'),
enanos: s('<circle cx="9" cy="13" r="3.5"/><circle cx="17.5" cy="7.5" r="1.5"/><circle cx="18" cy="16.5" r="1"/>'),
asteroides: s('<path d="M7 6l5-2 6 3 2 6-3 6-6 1-5-4-1-6z"/><circle cx="11" cy="10" r="1"/><circle cx="15" cy="14" r="1.2"/>'),
cometas: s('<circle cx="16.5" cy="7.5" r="2.5"/><path d="M14.5 9.5L4 20M13 7.5L5 13M17 11l-6 8"/>'),
constelaciones: s('<circle cx="5" cy="17" r="1"/><circle cx="10" cy="8" r="1"/><circle cx="16" cy="12" r="1"/><circle cx="20" cy="5" r="1"/><path d="M5.6 16.2l3.8-7.4M10.8 8.6l4.4 2.8M16.6 11.2l2.8-5.4"/>'),
naves: s('<rect x="9.5" y="9.5" width="5" height="5" rx="1"/><path d="M3 7.5l5 5M16 11.5l5 5" /><path d="M2.5 9.5l3-3 4 4-3 3zM14.5 13.5l3-3 4 4-3 3z"/><path d="M12 14.5v3M10 19.5a2.8 2.8 0 0 1 4 0"/>'),
capas: s('<path d="M12 3l9 5-9 5-9-5z"/><path d="M3 13l9 5 9-5"/>'),
config: s('<path d="M4 7h10M18 7h2M4 17h4M12 17h8"/><circle cx="16" cy="7" r="2"/><circle cx="10" cy="17" r="2"/>'),
search: s('<circle cx="11" cy="11" r="6.5"/><path d="M16 16l4.5 4.5"/>'),
close: s('<path d="M6 6l12 12M18 6L6 18"/>'),
menu: s('<path d="M4 7h16M4 12h16M4 17h16"/>'),
play: s('<path d="M8 5.5v13l10-6.5z" fill="currentColor"/>'),
pause: s('<path d="M8 5v14M16 5v14" stroke-width="2.4"/>'),
slower: s('<path d="M11 6l-6 6 6 6M19 6l-6 6 6 6"/>'),
faster: s('<path d="M13 6l6 6-6 6M5 6l6 6-6 6"/>'),
reset: s('<path d="M4.5 12a7.5 7.5 0 1 0 2.2-5.3M4.5 4.5v3.5H8"/>'),
travel: s('<path d="M5 19L19 5M9 5h10v10"/>'),
center: s('<circle cx="12" cy="12" r="3"/><path d="M12 2.5v4M12 17.5v4M2.5 12h4M17.5 12h4"/>'),
follow: s('<circle cx="12" cy="12" r="7.5" stroke-dasharray="3 3"/><circle cx="12" cy="12" r="2.5" fill="currentColor"/>'),
zoomin: s('<circle cx="11" cy="11" r="6.5"/><path d="M16 16l4.5 4.5M8.5 11h5M11 8.5v5"/>'),
zoomout: s('<circle cx="11" cy="11" r="6.5"/><path d="M16 16l4.5 4.5M8.5 11h5"/>'),
overview: s('<rect x="3.5" y="3.5" width="17" height="17" rx="2"/><circle cx="12" cy="12" r="1.5"/><ellipse cx="12" cy="12" rx="6" ry="3"/>'),
sun: s('<circle cx="12" cy="12" r="4"/><path d="M12 2.5v2.5M12 19v2.5M2.5 12H5M19 12h2.5M5.3 5.3l1.8 1.8M16.9 16.9l1.8 1.8M5.3 18.7l1.8-1.8M16.9 7.1l1.8-1.8"/>'),
map: s('<circle cx="12" cy="12" r="8.5"/><circle cx="12" cy="12" r="4.5"/><circle cx="12" cy="12" r="1" fill="currentColor"/><circle cx="16" cy="9.5" r="1.2" fill="currentColor"/>'),
calendar: s('<rect x="4" y="5.5" width="16" height="14.5" rx="2"/><path d="M4 10h16M8.5 3.5v4M15.5 3.5v4"/>'),
help: s('<circle cx="12" cy="12" r="9"/><path d="M9.6 9.5a2.5 2.5 0 1 1 3.4 2.3c-.7.3-1 .9-1 1.7M12 17h.01"/>'),
tour: s('<path d="M4 18c3-8 6 2 9-5s5-5 7-7"/><circle cx="20" cy="6" r="1.5"/>'),
};
})();
const NF = new Intl.NumberFormat('es-MX', { maximumFractionDigits: 2 });
const fmt = (x, d) => { const dd = d == null ? 2 : d; return new Intl.NumberFormat(I18N.loc(), { maximumFractionDigits: dd, minimumFractionDigits: 0 }).format(Math.abs(x) < 0.5 * Math.pow(10, -dd) ? 0 : x); };
const SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
function fmtSci(x, unit) {
if (x == null) return null;
const e = Math.floor(Math.log10(x)), m = x / Math.pow(10, e);
return fmt(m, 2) + ' × 10' + String(e).split('').map(c => SUP[c]).join('') + ' ' + unit;
}
function fmtKm(km) {
if (km < 1e6) return fmt(km, km < 100 ? 1 : 0) + ' km';
if (km < 1e9) return t('{0} millones de km', fmt(km / 1e6, km < 1e7 ? 2 : 1));
return t('{0} mil millones de km', fmt(km / 1e9, 2));
}
function fmtAU(au) { return fmt(au, au < 0.1 ? 4 : au < 10 ? 3 : 2) + ' ' + t('UA'); }
function fmtDur(days) {
if (days == null) return null;
const a = Math.abs(days);
if (a < 2) { const h = a * 24, hh = Math.floor(h), mm = Math.round((h - hh) * 60); return hh + ' h' + (mm ? ' ' + mm + ' min' : ''); }
if (a < 800) return t('{0} días', fmt(a, a < 30 ? 2 : 1));
return t('{0} años ({1} días)', fmt(a / 365.25, a / 365.25 < 100 ? 2 : 0), fmt(a, 0));
}
const norm = s => s.normalize('NFD').replace(/[\u0300-\u036f]/g, '').toLowerCase();
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const escNum = s => esc(s).replace(/^((?:[≈~]\s?)?[-−+]?\d[\d.,]*)/, '<span class="num">$1</span>');   // cifra inicial en tipografía técnica
const $ = s => document.querySelector(s);
const Info = {
orbitA(rb) { const o = rb.def.orbit; if (!o) return null; if (o.t === 'jpl') return o.el[0]; if (o.t === 'kep') return o.a; return null; },
periodDays(rb) {
const i = rb.def.info, o = rb.def.orbit;
if (i.year) return { v: i.year, k: 'ref' };
if (o && (o.t === 'moon' || o.t === 'luna')) return { v: o.P, k: 'ref' };
const a = this.orbitA(rb); if (a) return { v: 365.25 * Math.pow(a, 1.5), k: 'calc' };
return null;
},
fields(rb) {
const d = rb.def, i = d.info, o = d.orbit, out = { phys: [], orb: [], comp: [] };
const add = (sec, label, val, kind, note) => out[sec].push({ label, val, kind: kind || 'ref', note });
const isMoon = d.type === 'moon', par = rb.parent;
add('phys', 'Edad aproximada', i.age || null);
const diam = d.shape ? '≈ ' + fmt(2 * d.R, 1) + ' km (medio; forma irregular)' : fmt(2 * d.R, d.R < 50 ? 2 : 0) + ' km';
add('phys', 'Diámetro medio', d.R < 1 ? fmt(2 * d.R * 1000, 0) + ' m' + (d.shape ? ' (medio)' : '') : diam);
add('phys', 'Radio medio', d.R < 1 ? fmt(d.R * 1000, 0) + ' m' : fmt(d.R, d.R < 50 ? 2 : 1) + ' km');
add('phys', 'Masa', d.mass ? fmtSci(d.mass, 'kg') : null);
let g = i.grav, gk = 'ref';
if (g == null && d.mass) { g = G_CONST * d.mass / Math.pow(d.R * 1000, 2); gk = 'calc'; }
add('phys', 'Gravedad superficial', g != null ? fmt(g, g < 0.01 ? 5 : g < 1 ? 3 : 2) + ' m/s²' : null, gk);
if (d.mass) add('phys', 'Densidad media', fmt(d.mass * 1000 / (4 / 3 * Math.PI * Math.pow(d.R * 1e5, 3)), 2) + ' g/cm³', 'calc');   // masa en g / volumen en cm³
add('phys', 'Temperatura aproximada', i.temp || null);
let moons = i.moons;
if (typeof moons === 'number') moons = moons === 0 ? (isMoon ? 'Ninguna' : 'Ninguna conocida') : String(moons) + (moons > 20 ? ' (conocidas a marzo de 2026)' : '');
add('phys', 'Número de lunas', moons != null ? moons : null);
const P = this.periodDays(rb);
let rotTxt = i.rotTxt || null, rotK = 'ref';
if (!rotTxt && d.locked && P) { rotTxt = 'Síncrona: ' + fmtDur(P.v) + ' (igual a su periodo orbital)'; }
else if (!rotTxt && d.rotH) rotTxt = fmtDur(Math.abs(d.rotH) / 24) + (d.rotH < 0 ? ' (retrógrada)' : '');
let day = i.day || null, dayK = 'ref';
if (!day && d.locked && P) { day = '≈ ' + fmtDur(P.v) + ' (por su rotación síncrona)'; dayK = 'calc'; }
else if (!day && d.rotH && P && !isMoon) { const sd = 1 / Math.abs(1 / (d.rotH / 24) - Math.sign(d.rotH) / P.v); day = '≈ ' + fmtDur(sd); dayK = 'calc'; }
add('orb', 'Duración de un día', day, dayK);
add('orb', 'Periodo de rotación', rotTxt, rotK);
const yearLabel = isMoon ? 'Periodo orbital (alrededor de ' + par.def.name + ')' : d.type === 'star' ? 'Duración de un año' : 'Duración de un año (periodo orbital)';
add('orb', yearLabel, i.yearTxt || (P ? fmtDur(P.v) : null), i.yearTxt ? 'ref' : P ? P.k : 'ref');
let v = i.v, vk = 'ref';
if (v == null && !i.vTxt) {
if (isMoon && o.a && o.P) { v = TAU * o.a / (o.P * 86400); vk = 'calc'; }
else { const a = this.orbitA(rb); if (a && P) { v = TAU * a * AU_KM / (P.v * 86400); vk = 'calc'; } }
}
add('orb', 'Velocidad orbital media', i.vTxt || (v != null ? fmt(v, v < 10 ? 2 : 1) + ' km/s' : null), i.vTxt ? 'ref' : vk);
if (isMoon) add('orb', 'Distancia media a ' + par.def.name, fmtKm(o.a));
else if (d.type !== 'star') { const a = this.orbitA(rb); add('orb', 'Distancia media al Sol', a ? fmtKm(a * AU_KM) + ' (' + fmtAU(a) + ')' : null); }
if (o && (o.t === 'jpl' || o.t === 'kep')) {
const k = Astro.elementsAt(o, World.jd);
add('orb', 'Excentricidad', fmt(k.e, 4));
add('orb', 'Inclinación orbital', fmt(k.i, 2) + '° respecto a la eclíptica');
add('orb', 'Perihelio y afelio', fmtAU(k.a * (1 - k.e)) + ' y ' + fmtAU(k.a * (1 + k.e)), 'calc');
} else if (isMoon && o.t === 'moon') add('orb', 'Inclinación orbital', fmt(o.i, 2) + '° respecto al ecuador de ' + par.def.name + (o.i > 90 ? ' (retrógrada)' : ''));
else if (o && o.t === 'luna') add('orb', 'Inclinación orbital', '5.14° respecto a la eclíptica');
add('orb', 'Inclinación axial', i.tilt || null);
add('comp', 'Composición', i.comp || null);
add('comp', 'Tipo de atmósfera', i.atm || null);
return out;
},
craftFields(rb) {
const c = rb.def.craft, yrs = (World.jd - Astro.jdFromISO(c.launch)) / 365.25;
return [
{ label: 'Agencia', val: c.agency }, { label: 'Lanzamiento', val: c.launchTxt },
{ label: 'Tiempo en el espacio', val: yrs < 0 ? 'Aún no se había lanzado en la fecha simulada' : fmt(yrs, 1) + ' años', kind: 'calc' },
{ label: 'Ubicación', val: c.orbitTxt }, { label: 'Periodo orbital', val: c.periodTxt },
{ label: 'Dimensiones', val: c.dims }, { label: 'Masa', val: c.mass }, { label: 'Objetivo', val: c.goal },
{ label: 'Estado', val: t('{0} (según la información disponible a mediados de 2026)', c.status) },
];
},
craftNote(rb) {
return {
leo: 'Posición: órbita circular aproximada. La ubicación a lo largo de la órbita es ilustrativa porque no se usan datos de seguimiento en tiempo real.',
lpoint: 'Posición: órbita halo aproximada alrededor del punto de Lagrange; la amplitud y la fase son ilustrativas.',
drift: 'Posición: distancia extrapolada desde datos publicados por NASA/JPL, en una dirección fija aproximada.',
kep: 'Posición: órbita kepleriana aproximada; la orientación del plano orbital es una estimación.',
}[rb.def.orbit.t] + ' El tamaño de la nave en pantalla está muy ampliado.';
},
live(rb) {
const E = World.byId.tierra, out = [];
if (rb.isSun) {
const de = V.len(E.helio);
out.push({ label: 'Distancia a la Tierra', km: de * AU_KM, au: de });
return out;
}
const ds = V.len(rb.helio);
out.push({ label: 'Distancia actual al Sol', km: ds * AU_KM, au: ds });
if (rb !== E) { const de = V.dist(rb.helio, E.helio); out.push({ label: 'Distancia actual a la Tierra', km: de * AU_KM, au: de }); }
const o = rb.def.orbit;
if (o && (o.t === 'jpl' || o.t === 'kep')) {
const a = Astro.elementsAt(o, World.jd).a, v = Math.sqrt(1.32712440018e11 * (2 / (ds * AU_KM) - 1 / (a * AU_KM)));
out.push({ label: 'Velocidad orbital actual', txt: fmt(v, 2) + ' km/s' });
}
if (rb.isCraft && rb.def.orbit.t === 'leo') out.push({ label: 'Altitud sobre la Tierra', txt: fmtKm(rb.def.orbit.a - rb.parent.R) });
if (rb.isCraft && rb.def.orbit.t === 'drift') { const lt = V.dist(rb.helio, E.helio) * AU_KM / 299792.458 / 3600; out.push({ label: 'Tiempo de la señal desde la Tierra', txt: fmt(lt, 1) + ' h' }); }
if (rb.def.type === 'moon') { const dp = V.dist(rb.helio, rb.parent.helio) * AU_KM; out.push({ label: 'Distancia actual a ' + rb.parent.def.name, txt: fmtKm(dp) }); }
return out;
},
};
const UI = {
sel: null, hov: null, section: null, get reducedMotion() { return Settings.state.accessibility.reducedMotion; },
labels: {}, constLabels: [], tour: null, lastLive: 0, started: false, mmOn: true, toastT: 0,
init() {
this.buildRail(); this.buildCamTools(); this.buildTime(); this.buildModes(); this.buildSearch(); this.buildLabels(); this.bindKeys();
$('#btn-menu').innerHTML = ICON.menu;
$('#btn-menu').addEventListener('click', () => this.openSection(this.section ? null : 'explorar'));
$('#drawer-close').innerHTML = ICON.close; $('#drawer-close').addEventListener('click', () => this.openSection(null));
$('#info-close').innerHTML = ICON.close; $('#info-close').addEventListener('click', () => this.select(null));
$('#mm-toggle').innerHTML = ICON.map; $('#mm-toggle').addEventListener('click', () => this.toggleMinimap());
$('#help-btn').innerHTML = ICON.help; $('#help-btn').addEventListener('click', () => this.toggleHelp(true));
$('#help-close').addEventListener('click', () => this.toggleHelp(false));
$('#help-guide').addEventListener('click', () => { this.toggleHelp(false); this.toggleGuide(true); });
$('#guide-btn').addEventListener('click', () => this.toggleGuide(true));
$('#guide-x').addEventListener('click', () => this.toggleGuide(false));
$('#guide-ok').addEventListener('click', () => this.toggleGuide(false));
$('#guide').addEventListener('click', e => { if (e.target.id === 'guide') this.toggleGuide(false); });
addEventListener('keydown', e => { if (e.key === 'Escape' && $('#guide').classList.contains('open')) { e.stopImmediatePropagation(); this.toggleGuide(false); } }, true);
$('#help').addEventListener('click', e => { if (e.target.id === 'help') this.toggleHelp(false); });
$('#cine-exit').addEventListener('click', () => this.setMode('explore'));
$('#home-btn').addEventListener('click', () => App.home());
$('#fly-btn').addEventListener('click', () => Flight.openHangar());
this.initAudio();
this.buildTours();
this.initSettings();
this.initTourCine();
this.initAbout();
this.buildCredits();
this.initMusic();
$('#minimap').addEventListener('click', e => this.minimapClick(e));
this.tipEl = $('#tip');
document.addEventListener('mouseover', e => { const t = e.target.closest && e.target.closest('[data-tip]'); if (t) this.showTip(t); });
document.addEventListener('mouseout', e => { const t = e.target.closest && e.target.closest('[data-tip]'); if (t) this.tipEl.classList.remove('on'); });
document.addEventListener('focusin', e => { const t = e.target.closest && e.target.closest('[data-tip]'); if (t && t.matches(':focus-visible')) this.showTip(t); });
document.addEventListener('focusout', () => this.tipEl.classList.remove('on'));
document.addEventListener('pointerdown', () => this.tipEl.classList.remove('on'), true);
if (innerWidth < 760) { this.mmOn = false; }
this.refreshMinimap();
this.updateScaleBadge();
},
hint(msg, ms) { if (Settings.state.general.tips && !this.tourCine) this.toast(msg, ms); },
showTip(t) {
if (!Settings.state.general.tips) return;
const r = t.getBoundingClientRect(), el = this.tipEl;
el.textContent = t.dataset.tip; el.classList.add('on');
const w = el.offsetWidth, h = el.offsetHeight;
let x = r.left + r.width / 2 - w / 2, y = r.top - h - 8;
if (t.dataset.tipSide === 'right') { x = r.right + 10; y = r.top + r.height / 2 - h / 2; }
if (y < 6) y = r.bottom + 8;
el.style.transform = `translate(${clamp(x, 6, innerWidth - w - 6)}px, ${y}px)`;
},
toast(msg, ms) {
const t = $('#toast'); t.textContent = msg; t.classList.add('on');
clearTimeout(this.toastT); this.toastT = setTimeout(() => t.classList.remove('on'), ms || 3200);
},
userActed() { if (this.tour) this.stopTour(); },
SECTIONS: [
['explorar', 'Explorar'], ['sistema', 'Sistema Solar'], ['planetas', 'Planetas'], ['lunas', 'Lunas'], ['enanos', 'Planetas enanos'],
['asteroides', 'Asteroides'], ['cometas', 'Cometas'], ['naves', 'Naves espaciales'], ['constelaciones', 'Constelaciones'], ['capas', 'Capas'], ['config', 'Ajustes'],
],
buildRail() {
const r = $('#rail');
r.innerHTML = this.SECTIONS.map(([id, label]) => `<button class="rail-btn" data-sec="${id}" aria-label="${label}" data-tip="${label}" data-tip-side="right">${ICON[id]}<span>${label}</span></button>`).join('');
r.addEventListener('click', e => { const b = e.target.closest('[data-sec]'); if (!b) return; if (b.dataset.sec === 'config') { this.openSettings(); return; } this.openSection(this.section === b.dataset.sec ? null : b.dataset.sec); });
},
openSection(id) {
this.section = id;
document.querySelectorAll('.rail-btn').forEach(b => b.classList.toggle('on', b.dataset.sec === id));
const d = $('#drawer');
if (!id) { d.classList.remove('open'); d.setAttribute('aria-hidden', 'true'); return; }
$('#drawer-title').textContent = this.SECTIONS.find(s => s[0] === id)[1];
$('#drawer-body').innerHTML = (['planetas', 'enanos', 'asteroides', 'cometas', 'naves', 'lunas'].includes(id) ? this.only3DHTML() : '') + this.sectionHTML(id);
const o3 = $('#drawer-body [data-only3d]'); if (o3) o3.addEventListener('change', () => { this.only3D = o3.checked; this.applyOnly3D(); });
this.applyOnly3D();
d.classList.add('open'); d.setAttribute('aria-hidden', 'false');
this.bindSection(id);
if (id === 'config') this.refreshAudio();
if (innerWidth < 760) this.select(null, { keepCam: true });
},
badge(rb) { return `<span class="vb" data-vb="${rb.id}">${this.badgeInner(rb)}</span>`; },
badgeInner(rb) {
const v = Assets.visual(rb);
const cube = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.6l5.6 3.1v6.6L8 14.4l-5.6-3.1V4.7z M2.4 4.7L8 7.8l5.6-3.1 M8 7.8v6.6"/></svg>';
if (v.type === 'mesh') {
if (v.status === 'failed') return `<span class="b3d fail" data-tip="No se pudo cargar el modelo 3D importado; se muestra la representación estándar.">Modelo no disponible</span>`;
if (v.status === 'loading') return `<span class="b3d loading" data-tip="Descargando el modelo 3D importado…">${cube}3D · Cargando</span>`;
return `<span class="b3d" role="img" aria-label="Modelo 3D importado" data-tip="${v.status === 'loaded' ? 'Este objeto utiliza un modelo 3D importado.' : 'Este objeto utiliza un modelo 3D importado; se descarga al verlo de cerca.'}">${cube}3D</span>`;
}
if (v.type === 'texture') return `<span class="btex" data-tip="Superficie con textura importada sobre la esfera estándar (no es un modelo 3D).">Textura</span>`;
return '';
},
shipBadge(s, full) {
const st = PackedModels.state[s.id], has = !!MODEL_PACKS[s.id];
if (!has) return full ? 'Representación estándar (modelo procedural)' : '';
if (st === 'fail') return full ? '<span class="b3d fail">Modelo no disponible</span> Se usa la representación estándar.' : '<span class="b3d fail">Modelo no disponible</span>';
const cube = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.6l5.6 3.1v6.6L8 14.4l-5.6-3.1V4.7z M2.4 4.7L8 7.8l5.6-3.1 M8 7.8v6.6"/></svg>';
if (st === 'loading') return `<span class="b3d loading">${cube}3D · Cargando</span>` + (full ? ' Descargando el modelo 3D importado…' : '');
return `<span class="b3d" data-tip="Esta nave utiliza un modelo 3D importado.">${cube}3D</span>` + (full ? (st === 'ok' ? ' Modelo 3D importado, cargado correctamente.' : ' Modelo 3D importado; se descarga al abrir el hangar.') : '');
},
refreshBadges() {
document.querySelectorAll('[data-vbs]').forEach(el => { const s = SHIPS.find(x => x.id === el.dataset.vbs); if (s) { const h = this.shipBadge(s); if (el.innerHTML !== h) el.innerHTML = h; } });
document.querySelectorAll('[data-vbs-full]').forEach(el => { const s = SHIPS.find(x => x.id === el.dataset.vbsFull); if (s) { const h = this.shipBadge(s, true); if (el.innerHTML !== h) el.innerHTML = h; } });
document.querySelectorAll('[data-vb]').forEach(el => { const rb = World.byId[el.dataset.vb]; if (rb) { const h = this.badgeInner(rb); if (el.innerHTML !== h) el.innerHTML = h; } });
document.querySelectorAll('.row[data-go]').forEach(r => { const rb = World.byId[r.dataset.go]; if (rb) r.classList.toggle('has3d', Assets.visual(rb).type === 'mesh'); });
this.applyOnly3D();
},
only3D: false,
only3DHTML() { return `<div class="only3d-bar"><label class="mix-sw only3d"><input type="checkbox" data-only3d ${this.only3D ? 'checked' : ''}><span class="sw" aria-hidden="true"></span><span>Solo objetos con modelo 3D</span></label><span class="vb-legend"><span class="b3d">3D</span> modelo importado · <span class="btex">Textura</span> superficie importada</span></div>`; },
applyOnly3D() {
const d = $('#drawer-body'); if (!d) return;
d.classList.toggle('only3d-on', this.only3D);
d.querySelectorAll('.list').forEach(l => { const any = l.querySelector('.row.has3d'); l.classList.toggle('empty3d', this.only3D && !any); const h = l.previousElementSibling; if (h && h.tagName === 'H3') h.classList.toggle('empty3d', this.only3D && !any); });
let msg = d.querySelector('.only3d-none'); const none = this.only3D && !d.querySelector('.row.has3d');
if (none && !msg) { msg = document.createElement('p'); msg.className = 'note only3d-none'; msg.textContent = 'Ningún objeto de esta lista tiene todavía un modelo 3D importado.'; d.appendChild(msg); }
if (msg) msg.hidden = !none;
},
bodyRow(rb, extra) {
const c = rb.def.color;
const vt = Assets.visual(rb).type;
return `<button class="row${vt === 'mesh' ? ' has3d' : ''}" data-go="${rb.id}"><i class="dot" style="--c:${c}"></i><span class="row-main"><b>${esc(rb.def.short || rb.def.name)}${this.badge(rb)}</b><small>${esc(extra || TYPE_LABEL[rb.def.type])}</small></span><span class="row-go">${ICON.travel}</span></button>`;
},
sectionHTML(id) {
const W = World, by = t => W.rb.filter(r => r.def.type === t);
const distTxt = rb => fmtAU(V.len(rb.helio)) + ' del Sol ahora';
switch (id) {
case 'explorar': {
const feat = ['tierra', 'saturno', 'jupiter', 'marte', 'luna', 'jwst', 'pluton', 'halley'].map(i => W.byId[i]);
return `<p class="lead">Elige un destino o recorre el Sistema Solar con un viaje guiado. Arrastra para girar, usa la rueda para acercarte y haz clic en cualquier objeto.</p>
<button class="cta-row" id="tour-btn">${ICON.tour}<span><b>${this.tour ? 'Detener recorrido guiado' : 'Recorridos guiados'}</b><small>Planetas, lunas, exploración humana, asteroides y cometas</small></span></button>
<button class="cta-row" data-act="fly">${ICON.tour}<span><b>Pilotar una nave</b><small>Vuelo espacial libre con física de inercia</small></span></button>
<button class="cta-row" data-act="home">${ICON.overview}<span><b>Volver a la pantalla inicial</b><small>Ver de nuevo la presentación de SOLARIS</small></span></button>
<h3>Destinos destacados</h3><div class="list">${feat.map(rb => this.bodyRow(rb, rb.def.info.desc.split(':')[0].split('.')[0])).join('')}</div>
<h3>Controles</h3><dl class="keys"><dt>Arrastrar</dt><dd>Girar la cámara</dd><dt>Clic derecho o Mayús + arrastrar</dt><dd>Desplazar</dd><dt>Rueda o pellizco</dt><dd>Acercar y alejar</dd><dt>Doble clic</dt><dd>Viajar muy cerca</dd><dt>?</dt><dd>Todos los atajos</dd></dl>`;
}
case 'sistema': {
const sc = S.scale;
return `<p class="lead">8 planetas, 5 planetas enanos reconocidos, cientos de lunas y millones de asteroides y cometas giran alrededor del Sol. Aquí se representan ${W.rb.length} objetos individuales y más de 15 000 asteroides y objetos del cinturón de Kuiper.</p>
<button class="cta-row" data-act="overview">${ICON.overview}<span><b>Vista general</b><small>Encuadrar todo el sistema</small></span></button>
<h3>Escala</h3><div class="seg" role="radiogroup" aria-label="Escala">${Object.keys(SCALES).map(k => `<button role="radio" aria-checked="${k === sc}" class="${k === sc ? 'on' : ''}" data-scale="${k}">${SCALES[k].label}</button>`).join('')}</div>
<p class="note">${this.scaleNote(sc)}</p>
<h3>Qué es real y qué es aproximado</h3>
<ul class="facts"><li><b>Datos de referencia:</b> masas, radios, periodos y composición proceden de NASA y JPL.</li><li><b>Posiciones:</b> calculadas para la fecha simulada con elementos keplerianos aproximados de JPL (planetas, válidos entre 1800 y 2050) y elementos simplificados para el resto.</li><li><b>Modificaciones visuales:</b> en escala visual y educativa se amplían los tamaños y se comprimen las distancias. La iluminación se normaliza para que los planetas lejanos sean visibles.</li></ul>`;
}
case 'planetas': return `<div class="list">${by('planet').map(rb => this.bodyRow(rb, distTxt(rb))).join('')}</div>`;
case 'enanos': return `<h3>Planetas enanos</h3><div class="list">${by('dwarf').map(rb => this.bodyRow(rb, distTxt(rb))).join('')}</div><h3>Objetos transneptunianos</h3><div class="list">${by('tno').map(rb => this.bodyRow(rb, distTxt(rb))).join('')}</div>`;
case 'asteroides': return `<p class="lead">El cinturón principal está entre Marte y Júpiter. Observa los huecos de Kirkwood, despejados por resonancias con Júpiter, y los troyanos que comparten su órbita.</p><div class="list">${by('asteroid').map(rb => this.bodyRow(rb, distTxt(rb))).join('')}</div>`;
case 'cometas': return `<p class="lead">Las colas aparecen al acercarse al Sol: la de iones, azulada y recta, y la de polvo, curvada.</p><div class="list">${by('comet').map(rb => this.bodyRow(rb, distTxt(rb))).join('')}</div>
<label class="check"><input type="checkbox" data-layer="cometOrbits" ${S.layers.cometOrbits ? 'checked' : ''}><span>Mostrar trayectorias de cometas</span></label>`;
case 'naves': {
const g = (t, ids, note) => `<h3>${t}</h3><div class="list">${ids.map(id => { const rb = W.byId[id]; return this.bodyRow(rb, rb.hidden ? 'No existía en la fecha simulada' : rb.def.sub); }).join('')}</div>${note ? `<p class="note">${note}</p>` : ''}`;
return `<p class="lead">Telescopios, estaciones y sondas construidos por la humanidad. Su tamaño en pantalla está muy ampliado para que puedas verlos.</p>
<label class="check"><input type="checkbox" data-layer="craft" ${S.layers.craft ? 'checked' : ''}><span>Mostrar naves espaciales</span></label>` +
g('Órbita terrestre', ['iss', 'tiangong', 'hubble'], 'Órbitas circulares aproximadas: la posición exacta a lo largo de la órbita es ilustrativa.') +
g('Puntos de Lagrange', ['jwst', 'euclid', 'soho'], 'L1 y L2 están a unos 1.5 millones de km de la Tierra, hacia el Sol y en sentido opuesto.') +
g('Espacio profundo', ['parker', 'newhorizons', 'voyager2', 'voyager1'], 'Las sondas de escape se sitúan extrapolando su distancia y dirección desde datos publicados por NASA/JPL.');
}
case 'lunas': {
const groups = W.rb.filter(r => r.children.some(c => c.def.type === 'moon'));
return groups.map(p => `<h3>${esc(p.def.name)}</h3><div class="list">${p.children.filter(c => c.def.type === 'moon').map(m => this.bodyRow(m, 'a ' + fmtKm(m.def.orbit.a))).join('')}</div>`).join('');
}
case 'constelaciones': return `<label class="check"><input type="checkbox" data-layer="constLines" ${S.layers.constLines ? 'checked' : ''}><span>Líneas de constelaciones</span></label>
<label class="check"><input type="checkbox" data-layer="constNames" ${S.layers.constNames ? 'checked' : ''}><span>Nombres de constelaciones</span></label>
<p class="note">Las estrellas están a distancias enormes, por eso las constelaciones conservan su forma desde cualquier punto del Sistema Solar.</p>
<div class="list">${W.constellations.map(c => `<button class="row" data-const="${c.id}"><i class="dot star"></i><span class="row-main"><b>${esc(c.name)}</b><small>Mirar hacia esta constelación</small></span><span class="row-go">${ICON.travel}</span></button>`).join('')}</div>`;
case 'capas': {
const L = [['orbits', 'Órbitas planetarias'], ['moonOrbits', 'Órbitas de lunas'], ['cometOrbits', 'Trayectorias de cometas'], ['asteroids', 'Cinturón de asteroides'], ['kuiper', 'Cinturón de Kuiper'], ['moons', 'Lunas'],
['namesPlanets', 'Nombres de planetas'], ['namesMoons', 'Nombres de lunas'], ['namesSmall', 'Nombres de asteroides y cometas'], ['constNames', 'Nombres de constelaciones'], ['constLines', 'Líneas de constelaciones'],
['craft', 'Naves espaciales'], ['namesCraft', 'Nombres de naves espaciales'], ['orbitInfo', 'Información orbital en las etiquetas'], ['markers', 'Marcadores de objetos lejanos'], ['stars', 'Estrellas'], ['deepSky', 'Fondo profundo (Vía Láctea y nebulosas)']];
return L.map(([k, t]) => `<label class="check"><input type="checkbox" data-layer="${k}" ${S.layers[k] ? 'checked' : ''}><span>${t}</span></label>`).join('') +
`<label class="check"><input type="checkbox" id="chk-real" ${S.scale === 'real' ? 'checked' : ''}><span>Escala realista</span></label>`;
}
case 'config': return `<h3>Audio</h3>${mixerHTML()}<h3>Escala</h3><div class="seg" role="radiogroup" aria-label="Escala">${Object.keys(SCALES).map(k => `<button role="radio" aria-checked="${k === S.scale}" class="${k === S.scale ? 'on' : ''}" data-scale="${k}">${SCALES[k].label}</button>`).join('')}</div>
<h3>Calidad gráfica</h3><div class="seg" role="radiogroup" aria-label="Calidad">${['alta', 'media', 'baja'].map(k => `<button role="radio" aria-checked="${k === S.quality}" class="${k === S.quality ? 'on' : ''}" data-quality="${k}">${k[0].toUpperCase() + k.slice(1)}</button>`).join('')}</div>
<p class="note">Media y baja reducen la resolución interna, el antialiasing y el brillo difuso para mantener la fluidez.</p>
<h3>Luz en el lado nocturno</h3><input type="range" id="amb" min="0" max="0.25" step="0.005" value="${S.ambient}" aria-label="Luz ambiental">
<p class="note">Es una ayuda visual: en el espacio, el lado nocturno es casi negro.</p>
<label class="check"><input type="checkbox" id="chk-rotcap" ${S.rotCap ? 'checked' : ''}><span>Ralentizar la rotación a velocidades altas</span></label>
<p class="note">Evita el parpadeo cuando un día dura menos de un fotograma. La posición orbital sigue siendo exacta.</p>
<label class="check"><input type="checkbox" id="chk-mm" ${this.mmOn ? 'checked' : ''}><span>Mostrar minimapa</span></label>
<button class="cta-row" data-act="help">${ICON.help}<span><b>Atajos de teclado</b><small>Ver la lista completa</small></span></button>
<h3>Fuentes</h3><p class="note">NASA NSSDCA (fichas planetarias y de satélites), NASA Science, JPL Solar System Dynamics y Small-Body Database, elementos keplerianos aproximados de E. M. Standish (JPL). Recuento de lunas: Minor Planet Center y NASA, marzo de 2026.</p>
<h3>Modelos 3D</h3><ul class="facts">${[['jwst', 'Telescopio James Webb'], ['hubble', 'Telescopio Hubble'], ['iss', 'Estación Espacial Internacional'], ['voyager', 'Voyager 1 y 2']].map(([k, t]) => `<li><b>${t}:</b> ${this.creditHTML(MODEL_CREDITS[k])}</li>`).join('')}${Object.keys(PLANET_TEX).map(id => `<li><b>${BODY[id].name}:</b> ${this.creditHTML(BODY[id].credit)}</li>`).join('')}${SHIPS.filter(sh => MODEL_CREDITS[sh.id]).map(sh => `<li><b>Nave ${sh.name}:</b> ${this.creditHTML(MODEL_CREDITS[sh.id])}</li>`).join('')}</ul><p class="note">El resto de naves espaciales son modelos procedurales propios.</p>`;
}
return '';
},
bindSection(id) {
const b = $('#drawer-body');
b.querySelectorAll('[data-go]').forEach(el => el.addEventListener('click', () => { this.select(World.byId[el.dataset.go], { fly: true }); if (innerWidth < 760) this.openSection(null); }));
b.querySelectorAll('[data-layer]').forEach(el => el.addEventListener('change', () => this.setLayer(el.dataset.layer, el.checked)));
b.querySelectorAll('[data-scale]').forEach(el => el.addEventListener('click', () => this.setScale(el.dataset.scale)));
b.querySelectorAll('[data-quality]').forEach(el => el.addEventListener('click', () => { App.setQuality(el.dataset.quality); this.openSection(id); }));
b.querySelectorAll('[data-const]').forEach(el => el.addEventListener('click', () => {
const c = World.constellations.find(x => x.id === el.dataset.const);
this.setLayer('constLines', true); this.setLayer('constNames', true);
Cam.lookDir(c.dir); this.toast('Mirando hacia ' + c.name);
this.openSection('constelaciones');
}));
const act = b.querySelectorAll('[data-act]');
act.forEach(el => el.addEventListener('click', () => { if (el.dataset.act === 'overview') this.camAction('overview'); if (el.dataset.act === 'help') this.toggleHelp(true); if (el.dataset.act === 'home') App.home(); if (el.dataset.act === 'fly') Flight.openHangar(); }));
const tb = $('#tour-btn'); if (tb) tb.addEventListener('click', () => { if (this.tour) { this.stopTour(); this.openSection(innerWidth < 760 ? null : 'explorar'); } else this.toggleTours(true); });
const amb = $('#amb'); if (amb) amb.addEventListener('input', () => { S.ambient = +amb.value; });
const rc = $('#chk-rotcap'); if (rc) rc.addEventListener('change', () => { S.rotCap = rc.checked; });
const mm = $('#chk-mm'); if (mm) mm.addEventListener('change', () => this.toggleMinimap(mm.checked));
const cr = $('#chk-real'); if (cr) cr.addEventListener('change', () => this.setScale(cr.checked ? 'real' : 'visual'));
},
setLayer(k, v) {
S.layers[k] = v;
const L = Settings.state.simulation.layers; if (k in L && L[k] !== v) { L[k] = v; Settings.save(); }
document.querySelectorAll(`[data-layer="${k}"]`).forEach(el => { el.checked = v; });
if (k === 'moons' && !v && this.sel && this.sel.def.type === 'moon') this.select(this.sel.parent, { fly: true });
},
scaleNote(k) {
return {
visual: 'Los tamaños de los cuerpos están ampliados y las distancias comprimidas para que todo quepa en pantalla. Las posiciones relativas y el orden son correctos.',
real: 'Distancias y tamaños proporcionales entre sí. Los planetas se vuelven diminutos: usa las etiquetas, los marcadores o el buscador para encontrarlos.',
edu: 'Separación casi uniforme entre órbitas y tamaños muy ampliados. Útil para comparar los cuerpos de un vistazo, no para medir distancias.',
}[k];
},
setScale(k) {
if (k === S.scale) return;
const oldR = Cam.focus ? Cam.focus.rS : ScaleState.get('overview');
App.scaleRatio = Cam.dist / oldR;
S.scale = k; ScaleState.set(k, true);
if (Settings.state.simulation.scale !== k) { Settings.state.simulation.scale = k; Settings.save(); }
document.querySelectorAll('[data-scale]').forEach(el => { el.classList.toggle('on', el.dataset.scale === k); el.setAttribute('aria-checked', el.dataset.scale === k); });
const cr = $('#chk-real'); if (cr) cr.checked = k === 'real';
const note = document.querySelector('#drawer-body .note'); if (note && this.section === 'sistema') note.textContent = this.scaleNote(k);
this.updateScaleBadge();
this.toast(SCALES[k].badge);
},
updateScaleBadge() { $('#scalebadge').textContent = SCALES[S.scale].badge; },
CAM: [['explore', 'travel', 'Explorar el objeto seleccionado'], ['center', 'center', 'Centrar'], ['follow', 'follow', 'Seguir'], ['in', 'zoomin', 'Acercar'], ['out', 'zoomout', 'Alejar'], ['overview', 'overview', 'Vista general'], ['sun', 'sun', 'Volver al Sol']],
buildCamTools() {
const c = $('#camtools');
c.innerHTML = this.CAM.map(([a, ic, t]) => `<button class="tool" data-cam="${a}" aria-label="${t}" data-tip="${t}">${ICON[ic]}</button>`).join('');
c.addEventListener('click', e => { const b = e.target.closest('[data-cam]'); if (b) this.camAction(b.dataset.cam); });
this.refreshCamButtons();
},
camAction(a) {
this.userActed();
const s = this.sel || Cam.focus;
if (a === 'explore') { if (s) Cam.travel(s); else this.toast(matchMedia('(pointer: coarse)').matches ? 'Selecciona primero un objeto: tócalo o búscalo.' : 'Selecciona primero un objeto: haz clic sobre él o búscalo.'); }
if (a === 'center') { if (s) { if (Cam.focus !== s) Cam.travel(s, { keepAngles: true, dist: clamp(Cam.dist, s.rS * 1.3, Cam.maxDist()) }); else Cam.center(); } else this.toast('Selecciona un objeto para centrarlo.'); }
if (a === 'follow') { if (s) { if (Cam.focus === s && Cam.follow) { Cam.follow = false; this.toast('Seguimiento desactivado'); } else { if (Cam.focus !== s) Cam.travel(s, { keepAngles: true, dist: clamp(Cam.dist, s.rS * 1.3, Cam.maxDist()) }); Cam.follow = true; this.toast('Siguiendo a ' + s.def.name); } } else this.toast('Selecciona un objeto para seguirlo.'); }
if (a === 'in') Cam.zoom(0.55);
if (a === 'out') Cam.zoom(1.8);
if (a === 'overview') { Cam.travelTo([0, 0, 0], ScaleState.get('overview'), 0.42); }
if (a === 'sun') this.select(World.sun, { fly: true });
this.refreshCamButtons();
},
refreshCamButtons() {
const f = document.querySelector('[data-cam="follow"]');
if (f) { const on = !!(Cam.focus && Cam.follow && (!this.sel || Cam.focus === this.sel)); f.classList.toggle('on', on); f.setAttribute('aria-pressed', on); }
const ip = $('#info-follow'); if (ip) { const on = !!(this.sel && Cam.focus === this.sel && Cam.follow); ip.classList.toggle('on', on); ip.querySelector('span').textContent = on ? 'Siguiendo' : 'Seguir'; }
},
buildTime() {
$('#t-slower').innerHTML = ICON.slower; $('#t-faster').innerHTML = ICON.faster; $('#t-reset').innerHTML = ICON.reset;
$('#t-slower').addEventListener('click', () => this.speed(-1));
$('#t-faster').addEventListener('click', () => this.speed(1));
$('#t-play').addEventListener('click', () => this.togglePause());
$('#t-reset').addEventListener('click', () => { Time.goLive(); this.refreshTime(); this.toast('Fecha y hora actuales en tiempo real'); });
$('#t-now').addEventListener('click', () => { Time.goLive(); this.refreshTime(); this.toast('Fecha y hora actuales en tiempo real'); });
const pick = $('#t-pick');
pick.addEventListener('change', () => {
if (!pick.value) return;
this.setSimDate(Astro.jdFromISO(pick.value) + 0.5);          // único punto de cambio de fecha
});
$('#t-cal').innerHTML = ICON.calendar;
$('#t-cal').addEventListener('click', () => { try { pick.showPicker(); } catch (e) { pick.focus(); } });
const sp = $('#t-speed'), menu = $('#speed-menu');
menu.innerHTML = Time.speeds.map((s, i) => `<button role="menuitemradio" data-i="${i}">${s === 1 ? 'Tiempo real' : '×' + s.toLocaleString(I18N.loc()).replace(/,/g, ' ')}</button>`).join('');
sp.addEventListener('click', () => { const open = !menu.classList.contains('open'); menu.classList.toggle('open', open); sp.setAttribute('aria-expanded', open); });
menu.addEventListener('click', e => { const b = e.target.closest('[data-i]'); if (!b) return; Time.idx = +b.dataset.i; Time.paused = false; Time.live = false; menu.classList.remove('open'); this.refreshTime(); });
document.addEventListener('click', e => { if (!e.target.closest('#t-speed') && !e.target.closest('#speed-menu')) menu.classList.remove('open'); });
this.refreshTime();
},
speed(d) { Time.idx = clamp(Time.idx + d, 0, Time.speeds.length - 1); Time.paused = false; Time.live = false; this.refreshTime(); },
togglePause() { Time.paused = !Time.paused; Time.live = false; this.refreshTime(); },
dateStr(jd) { return Astro.dateFromJD(jd).toLocaleDateString(I18N.loc(), { day: 'numeric', month: 'short', year: 'numeric' }); },
refreshTime() {
const p = $('#t-play');
p.innerHTML = Time.paused ? ICON.play : ICON.pause;
p.setAttribute('aria-label', Time.paused ? 'Reanudar' : 'Pausar'); p.dataset.tip = Time.paused ? 'Reanudar (espacio)' : 'Pausar (espacio)';
$('#t-speed').textContent = Time.paused ? 'En pausa' : Time.label();
$('#t-speed').dataset.tip = Time.paused ? 'Simulación detenida' : Time.rateText();
document.querySelectorAll('#speed-menu [data-i]').forEach(b => b.setAttribute('aria-checked', +b.dataset.i === Time.idx));
},
tickTime() {
const d = Astro.dateFromJD(Time.jd);
if (isNaN(d.getTime())) return;
$('#t-date').textContent = d.toLocaleDateString(I18N.loc(), { day: 'numeric', month: 'short', year: 'numeric' });
const nowT = performance.now(); if (!this._tlT || nowT - this._tlT > 250) { this._tlT = nowT; TL.refresh(); }   // indicador de simulación (4 veces por segundo)
const live = Time.live && !Time.paused && Time.idx === 0;
$('#t-clock').textContent = d.toLocaleTimeString(I18N.loc(), { hour: '2-digit', minute: '2-digit', hour12: false }) + (live ? ' · en vivo' : ' hora local');
$('#t-clock').dataset.tip = d.toLocaleString(I18N.loc(), { day: 'numeric', month: 'short', year: 'numeric', hour: '2-digit', minute: '2-digit', timeZone: 'UTC', hour12: false }) + ' UTC';
const capped = World.rb.some(r => r.capped && r.proj.rpx > 20);
$('#rot-flag').classList.toggle('on', !!capped);
},
MODES: [['explore', 'Exploración', 'Navegación libre'], ['follow', 'Seguimiento', 'La cámara acompaña al objeto seleccionado'], ['orbits', 'Órbitas', 'Resalta todas las trayectorias'], ['edu', 'Educativo', 'Muestra nombres y datos orbitales'], ['cine', 'Cine', 'Oculta los controles y gira lentamente']],
buildModes() {
const m = $('#modes');
m.innerHTML = this.MODES.map(([k, t, tip]) => `<button role="radio" data-mode="${k}" aria-checked="${k === S.mode}" class="${k === S.mode ? 'on' : ''}" data-tip="${tip}">${t}</button>`).join('');
m.addEventListener('click', e => { const b = e.target.closest('[data-mode]'); if (b) this.setMode(b.dataset.mode); });
const fade = () => { const max = m.scrollWidth - m.clientWidth; m.classList.toggle('fade-l', max > 2 && m.scrollLeft > 2); m.classList.toggle('fade-r', max > 2 && m.scrollLeft < max - 2); };
m.addEventListener('scroll', fade, { passive: true }); addEventListener('resize', fade); requestAnimationFrame(fade);
},
setMode(k) {
const prev = S.mode; S.mode = k;
document.querySelectorAll('[data-mode]').forEach(b => { b.classList.toggle('on', b.dataset.mode === k); b.setAttribute('aria-checked', b.dataset.mode === k); });
{ const m = $('#modes'), on = m.querySelector('.on'); if (on && m.scrollWidth > m.clientWidth) m.scrollTo({ left: on.offsetLeft - (m.clientWidth - on.offsetWidth) / 2, behavior: 'smooth' }); }
const app = $('#app');
app.classList.toggle('cine', k === 'cine');
Cam.auto = k === 'cine' ? 0.035 : 0;
if (prev === 'edu' && k !== 'edu' && this.eduSaved) { Object.assign(S.layers, this.eduSaved); this.eduSaved = null; }
if (k === 'edu') { this.eduSaved = { namesMoons: S.layers.namesMoons, namesSmall: S.layers.namesSmall, orbitInfo: S.layers.orbitInfo }; Object.assign(S.layers, { namesMoons: true, namesSmall: true, orbitInfo: true }); this.toast('Modo educativo: nombres y datos orbitales visibles'); }
if (k === 'follow') { const t = this.sel || World.byId.tierra; this.select(t, { fly: Cam.focus !== t }); Cam.follow = true; this.toast('Siguiendo a ' + t.def.name); }
if (k === 'orbits') { this.toast('Órbitas resaltadas'); if (!this.sel) Cam.travelTo([0, 0, 0], ScaleState.get('overview') * 0.85, 1.1); }
if (k === 'cine') { this.select(null, { keepCam: true }); this.openSection(null); this.hint(I18N.t('Modo cine: pulsa Esc para salir'), 2400); }
if (this.section) this.openSection(this.section);
this.refreshCamButtons();
},
buildSearch() {
$('#search-ico').innerHTML = ICON.search;
const q = $('#q'), res = $('#results');
this.searchIdx = -1;
const run = () => {
const t = norm(q.value.trim());
if (!t) { res.classList.remove('open'); res.innerHTML = ''; q.setAttribute('aria-expanded', 'false'); return; }
const hits = World.rb.map(rb => {
const names = [rb.def.name, rb.def.short || '', ...(rb.def.aka || [])].map(norm);
let sc = 99;
names.forEach(n => { if (n === t) sc = Math.min(sc, 0); else if (n.startsWith(t)) sc = Math.min(sc, 1); else if (n.includes(t)) sc = Math.min(sc, 2); });
return { rb, sc };
}).filter(h => h.sc < 99).sort((a, b) => a.sc - b.sc || a.rb.def.name.localeCompare(b.rb.def.name)).slice(0, 8);
this.hits = hits.map(h => h.rb); this.searchIdx = hits.length ? 0 : -1;
res.innerHTML = hits.length ? hits.map((h, i) => `<li role="option" id="opt-${i}" aria-selected="${i === 0}" data-go="${h.rb.id}"><i class="dot" style="--c:${h.rb.def.color}"></i><b>${esc(h.rb.def.name)}</b>${this.badge(h.rb)}<small>${esc(TYPE_LABEL[h.rb.def.type])}${h.rb.parent && h.rb.def.type === 'moon' ? ' de ' + esc(h.rb.parent.def.name) : ''}</small></li>`).join('')
: `<li class="empty">Sin resultados para "${esc(q.value)}". Prueba con un planeta, luna, asteroide o cometa.</li>`;
res.classList.add('open'); q.setAttribute('aria-expanded', 'true');
};
const go = rb => { if (!rb) return; this.select(rb, { fly: true }); q.value = ''; run(); q.blur(); };
q.addEventListener('input', run);
q.addEventListener('keydown', e => {
if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
e.preventDefault(); if (!this.hits || !this.hits.length) return;
this.searchIdx = (this.searchIdx + (e.key === 'ArrowDown' ? 1 : -1) + this.hits.length) % this.hits.length;
res.querySelectorAll('li').forEach((li, i) => li.setAttribute('aria-selected', i === this.searchIdx));
q.setAttribute('aria-activedescendant', 'opt-' + this.searchIdx);
} else if (e.key === 'Enter') { e.preventDefault(); if (this.hits && this.searchIdx >= 0) go(this.hits[this.searchIdx]); }
else if (e.key === 'Escape') { q.value = ''; run(); q.blur(); }
});
res.addEventListener('mousedown', e => { const li = e.target.closest('[data-go]'); if (li) { e.preventDefault(); go(World.byId[li.dataset.go]); } });
q.addEventListener('blur', () => setTimeout(() => res.classList.remove('open'), 120));
q.addEventListener('focus', () => { if (q.value) run(); });
},
select(rb, opts) {
opts = opts || {};
if (this.tour && !opts.tour) this.stopTour(true);
const prev = this.sel;
this.sel = rb;
const panel = $('#info');
if (!rb) { panel.classList.remove('open'); panel.setAttribute('aria-hidden', 'true'); this.refreshCamButtons(); return; }
if (rb.def.type === 'moon' && !S.layers.moons) this.setLayer('moons', true);
if (rb.isCraft && !S.layers.craft) this.setLayer('craft', true);
if (rb.hidden) { this.sel = prev; this.toast(rb.def.name + ' no estaba en esta posición en la fecha simulada (lanzamiento: ' + rb.def.craft.launchTxt.split(' (')[0] + ').', 4200); return; }
if (prev !== rb) SFX.play('select');
if (opts.fly) Cam.travel(rb, { close: opts.close });
if (prev !== rb || !panel.classList.contains('open')) this.renderInfo(rb);
if (S.mode !== 'cine') { panel.classList.add('open'); panel.setAttribute('aria-hidden', 'false'); }
if (innerWidth < 760) this.openSection(null);
this.refreshCamButtons();
},
setHover(rb) {
if (this.hov === rb) return;
if (this.hov) this.hov.hover = 0;
this.hov = rb; if (rb) { rb.hover = 1; SFX.hover(); }
$('#gl').style.cursor = rb ? 'pointer' : 'grab';
},
renderInfo(rb) {
const d = rb.def, P = $('#info-body'), basic = Settings.state.interface.info === 'basic';
let F = rb.isCraft ? null : Info.fields(rb);
if (F && basic) {
const keep = l => /^(Diámetro medio|Gravedad superficial|Temperatura aproximada|Número de lunas|Duración de un día|Duración de un año|Periodo orbital|Distancia media)/.test(l);
F = { phys: F.phys.filter(f => keep(f.label)), orb: F.orb.filter(f => keep(f.label)), comp: [] };
}
const craftF = rb.isCraft ? Info.craftFields(rb).filter(f => !basic || /^(Agencia|Lanzamiento|Ubicación|Objetivo|Estado)/.test(f.label)) : null;
const tag = k => k === 'calc' ? '<em class="tag calc" data-tip="Calculado a partir de otros datos de referencia">calc.</em>' : k === 'sim' ? '<em class="tag sim" data-tip="Aproximación de la simulación">sim.</em>' : '';
const row = f => `<div class="kv"><dt>${esc(f.label)}</dt><dd>${f.val == null ? '<span class="na">No disponible</span>' : escNum(f.val)} ${f.val == null ? '' : tag(f.kind)}</dd></div>`;
let sub = TYPE_LABEL[d.type];
if (d.sub) sub = d.sub;
if (d.type === 'moon') sub = 'Luna de ' + rb.parent.def.name;
if (rb.isCraft) sub = d.sub;
const moons = rb.children.filter(c => c.def.type === 'moon');
P.innerHTML = `
<header class="info-head">
<span class="kind"><i class="dot" style="--c:${d.color}"></i>${esc(sub)}</span>
<h2 id="info-title">${esc(d.name)}</h2>
<p class="desc">${esc(d.info.desc)}</p>
<div class="info-actions">
<button class="pill" id="info-travel">${ICON.travel}<span>Viajar</span></button>
<button class="pill" id="info-follow">${ICON.follow}<span>Seguir</span></button>
<button class="pill" id="info-close2">${ICON.zoomin}<span>Acercar</span></button>
${Compare.pillHTML(rb)}
<button class="pill" id="info-plan" data-plan="${rb.id}">Planificar viaje</button>
</div>
${rb.isCraft && rb.def.orbit.t === 'lpoint' ? `<p class="rel">En el punto L${rb.def.orbit.L} Sol-Tierra, cerca de <button class="link" data-go="${rb.parent.id}" style="--c:${rb.parent.def.color || '#94a0b6'}">${esc(rb.parent.def.name)}</button></p>` : ''}
${rb.parent && !rb.parent.isSun && !(rb.isCraft && rb.def.orbit.t === 'lpoint') ? `<p class="rel">Orbita a <button class="link" data-go="${rb.parent.id}" style="--c:${rb.parent.def.color || '#94a0b6'}">${esc(rb.parent.def.name)}</button></p>` : ''}
${moons.length ? `<p class="rel">Lunas incluidas: ${moons.map(m => `<button class="link" data-go="${m.id}" style="--c:${m.def.color || '#94a0b6'}">${esc(m.def.name)}</button>`).join('')}</p>` : ''}
</header>
<section class="live"><h3>Ahora mismo <em class="tag sim" data-tip="Calculado por la simulación para la fecha mostrada; no es una efeméride de alta precisión">sim.</em></h3><div id="live-rows"></div></section>
${rb.isCraft ? `<section><h3>Misión</h3><dl>${craftF.map(row).join('')}</dl></section><p class="note">${Info.craftNote(rb)}</p>` : `<section><h3>Datos físicos</h3><dl>${F.phys.map(row).join('')}</dl></section>
<section><h3>Rotación y órbita</h3><dl>${F.orb.map(row).join('')}</dl></section>
${F.comp.length ? `<section><h3>Composición</h3><dl>${F.comp.map(row).join('')}</dl></section>` : ''}`}
${basic ? '<p class="note">Vista básica: activa la información científica avanzada en Ajustes › Interfaz para ver todos los datos.</p>' : ''}
${d.info.feats && d.info.feats.length ? `<section><h3>Características principales</h3><ul class="facts">${d.info.feats.map(x => `<li>${esc(x)}</li>`).join('')}</ul></section>` : ''}
${d.info.facts && d.info.facts.length ? `<section><h3>Curiosidades</h3><ul class="facts">${d.info.facts.map(x => `<li>${esc(x)}</li>`).join('')}</ul></section>` : ''}
<footer class="src"><p><b>Fuente:</b> ${esc(d.src)}</p>${d.credit && World.hd[d.model || d.id] ? `<p><b>Modelo 3D:</b> ${UI.creditHTML(d.credit)}</p>` : ''}<p class="legend"><em class="tag calc">calc.</em> derivado de otros datos <em class="tag sim">sim.</em> aproximación de la simulación</p>${this.scaleWarn(rb)}</footer>`;
P.scrollTop = 0;
$('#info-travel').addEventListener('click', () => { this.userActed(); Cam.travel(rb); });
$('#info-follow').addEventListener('click', () => { this.userActed(); if (Cam.focus === rb && Cam.follow) Cam.follow = false; else { if (Cam.focus !== rb) Cam.travel(rb, { keepAngles: true, dist: clamp(Cam.dist, rb.rS * 1.3, Cam.maxDist()) }); Cam.follow = true; } this.refreshCamButtons(); });
$('#info-close2').addEventListener('click', () => { this.userActed(); Cam.travel(rb, { close: true }); });
P.querySelectorAll('[data-go]').forEach(el => el.addEventListener('click', () => this.select(World.byId[el.dataset.go], { fly: true })));
this.updateLive(true);
this.refreshCamButtons();
},
creditHTML(c) {
if (!c.url) return esc(c.note.charAt(0).toUpperCase() + c.note.slice(1)) + '.';
return `<a href="${c.url}" target="_blank" rel="noopener">"${esc(c.title)}"</a> de ${esc(c.author)} (Sketchfab), licencia <a href="${c.licenseUrl}" target="_blank" rel="noopener">${esc(c.license)}</a>; ${esc(c.note)}.`;
},
scaleWarn(rb) {
if (S.scale === 'real') return '';
return `<p class="legend">En la escala ${SCALES[S.scale].label.toLowerCase()} el tamaño de ${esc(rb.def.name)} y sus distancias en pantalla están modificados; las cifras de esta ficha son reales.</p>`;
},
updateLive(force) {
const now = performance.now();
if (!this.sel || (!force && now - this.lastLive < 250)) return;
this.lastLive = now;
const rows = Info.live(this.sel);
const el = $('#live-rows'); if (!el) return;
el.innerHTML = rows.map(r => `<div class="live-row"><span>${esc(r.label)}</span><b>${escNum(r.txt || fmtKm(r.km))}</b>${r.au != null && r.au >= 0.001 ? `<small>${escNum(fmtAU(r.au))}</small>` : ''}</div>`).join('');
},
buildLabels() {
const host = $('#labels');
const all = World.systems ? [].concat(...Object.values(World.systems).map(s => s.rb)) : World.rb;   // etiquetas de todos los sistemas
all.forEach(rb => {
const el = document.createElement('div');
el.className = 'lbl lbl-' + rb.def.type;
el.innerHTML = `<span class="n">${esc(rb.def.short || rb.def.name)}</span><span class="o"></span>`;
host.appendChild(el); this.labels[rb.id] = { el, n: el.querySelector('.n'), o: el.querySelector('.o'), shown: false, w: 0 };
});
World.constellations.forEach(c => {
const el = document.createElement('div'); el.className = 'clbl'; el.textContent = c.name; host.appendChild(el);
this.constLabels.push({ c, el });
});
},
labelAllowed(rb, forPick) {
const L = S.layers, t = rb.def.type;
if (t === 'moon' && !L.moons) return false;
if (t === 'craft' && (!L.craft || rb.hidden)) return false;
if (forPick) return true;
if (S.mode === 'cine') return false;          // modo cine: ninguna etiqueta
if (rb === this.sel || rb === this.hov) return true;
if (t === 'star' || t === 'planet' || t === 'blackhole') return L.namesPlanets;
if (t === 'moon') return L.namesMoons;
if (t === 'dwarf') return L.namesPlanets || L.namesSmall;
if (t === 'craft') return L.namesCraft;
return L.namesSmall;
},
updateLabels() {
Sci.labels(); Planner.marks(); Systems.frame(); Missions.frame();   // capas científicas, ruta del planificador y avisos del sistema activo
const placed = [], W = innerWidth, H = innerHeight;
const pri = rb => (rb === this.sel ? 0 : rb === this.hov ? 1 : rb.isSun ? 2 : rb.def.type === 'planet' ? 3 : rb.def.type === 'dwarf' ? 4 : rb.def.type === 'moon' ? 6 : rb.isCraft ? 7 : 5);
const list = World.rb.slice().sort((a, b) => pri(a) - pri(b));
const showInfo = S.layers.orbitInfo || S.mode === 'edu';
for (const rb of list) {
const L = this.labels[rb.id], p = rb.proj;
let ok = p.on && this.labelAllowed(rb) && p.x > -60 && p.x < W + 60 && p.y > -30 && p.y < H + 30;
if (ok && (rb.def.type === 'moon' || (rb.isCraft && !rb.parent.isSun)) && rb !== this.sel && rb !== this.hov) {
const pp = rb.parent.proj, sep = Math.hypot(p.x - pp.x, p.y - pp.y);
ok = pp.on && sep > Math.max(26, pp.rpx + 14);
}
if (ok && rb !== this.sel) { // oculto tras un cuerpo más cercano
for (const o of World.rb) { if (o === rb || !o.proj.on || o.proj.rpx < 4) continue; if (o.proj.w < p.w && Math.hypot(o.proj.x - p.x, o.proj.y - p.y) < o.proj.rpx * 0.98) { ok = false; break; } }
}
const off = Math.max(p.rpx, 3) * 0.72 + 8;
const x = p.x + off, y = p.y - off * 0.62;
if (ok) {
if (!L.w) L.w = L.el.offsetWidth || 60;
const box = [x, y - 9, x + L.w, y + (showInfo ? 22 : 9)];
if (rb !== this.sel && placed.some(b => !(box[2] < b[0] || box[0] > b[2] || box[3] < b[1] || box[1] > b[3]))) ok = false;
else placed.push(box);
}
if (ok) {
L.el.style.transform = `translate(${x.toFixed(1)}px, ${(y - 9).toFixed(1)}px)`;
if (showInfo) { const txt = rb.isSun ? '' : (rb.def.type === 'moon' || (rb.isCraft && !rb.parent.isSun)) ? fmtKm(V.dist(rb.helio, rb.parent.helio) * AU_KM) : fmtAU(V.len(rb.helio)); if (L.o.textContent !== txt) L.o.textContent = txt; }
else if (L.o.textContent) L.o.textContent = '';
}
if (ok !== L.shown) { L.el.classList.toggle('on', ok); L.shown = ok; }
L.el.classList.toggle('sel', rb === this.sel); L.el.classList.toggle('hov', rb === this.hov);
}
const cn = S.layers.constNames && S.mode !== 'cine';
for (const c of this.constLabels) {
let ok = false;
if (cn) {
const v = World.view, d = c.c.dir;
const vx = v[0] * d[0] + v[4] * d[1] + v[8] * d[2], vy = v[1] * d[0] + v[5] * d[1] + v[9] * d[2], vz = v[2] * d[0] + v[6] * d[1] + v[10] * d[2];
if (vz < -0.05) {
const P = World.proj, x = ((P[0] * vx / -vz - P[8]) * 0.5 + 0.5) * W, y = (1 - ((P[5] * vy / -vz - P[9]) * 0.5 + 0.5)) * H;
ok = x > 0 && x < W && y > 0 && y < H;
if (ok) c.el.style.transform = `translate(${x.toFixed(0)}px, ${y.toFixed(0)}px) translate(-50%, -50%)`;
}
}
c.el.classList.toggle('on', ok);
}
const ret = (el, rb) => {
if (!rb || !rb.proj.on) { el.classList.remove('on'); return; }
const s = Math.max(rb.proj.rpx * 2 + 16, 26);
el.style.transform = `translate(${(rb.proj.x - s / 2).toFixed(1)}px, ${(rb.proj.y - s / 2).toFixed(1)}px)`;
el.style.width = el.style.height = s + 'px';
el.classList.add('on');
};
ret($('#reticle'), S.mode === 'cine' ? null : this.sel);
ret($('#reticle-hover'), this.hov !== this.sel ? this.hov : null);
},
toggleMinimap(v) {
this.mmOn = v == null ? !this.mmOn : v;
if (Settings.state.interface.minimap !== this.mmOn) { Settings.state.interface.minimap = this.mmOn; Settings.save(); }
$('#minimap-wrap').classList.toggle('off', !this.mmOn);
$('#app').classList.toggle('mm-on', this.mmOn);
$('#mm-toggle').classList.toggle('on', this.mmOn);
$('#mm-toggle').setAttribute('aria-pressed', this.mmOn);
const c = $('#chk-mm'); if (c) c.checked = this.mmOn;
},
refreshMinimap() { this.toggleMinimap(this.mmOn); },
mmMap(x, z, R) { const r = Math.hypot(x, z); if (r < 1e-9) return [0, 0]; const k = R * Math.log(1 + r / 0.35) / Math.log(1 + 45 / 0.35) / r; return [x * k, z * k]; },
drawMinimap() {
if (!this.mmOn) return;
const cv = $('#minimap'), dpr = Math.min(devicePixelRatio || 1, 2), S2 = cv.clientWidth;
if (!S2) return;
if (cv.width !== S2 * dpr) { cv.width = cv.height = S2 * dpr; }
const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
ctx.clearRect(0, 0, S2, S2);
const c = S2 / 2, R = S2 / 2 - 8;
ctx.save(); ctx.translate(c, c);
ctx.lineWidth = 1;
this.mmPts = [];
for (const rb of World.rb) {
const t = rb.def.type; if (!(t === 'planet' || t === 'dwarf') || !rb.orbitPts) continue;
ctx.strokeStyle = rb === this.sel ? 'rgba(242,200,121,0.9)' : t === 'planet' ? 'rgba(150,180,230,0.28)' : 'rgba(194,174,140,0.18)';
ctx.beginPath();
for (let i = 0; i < rb.orbitPts.length; i += 6) { const p = rb.orbitPts[i], m = this.mmMap(p[0], p[2], R); if (i) ctx.lineTo(m[0], m[1]); else ctx.moveTo(m[0], m[1]); }
ctx.closePath(); ctx.stroke();
}
ctx.fillStyle = '#ffd27a'; ctx.beginPath(); ctx.arc(0, 0, 3, 0, TAU); ctx.fill();
for (const rb of World.rb) {
const t = rb.def.type; if (!(t === 'planet' || t === 'dwarf')) continue;
const h = Astro.eclToScene(rb.helio), m = this.mmMap(h[0], h[2], R);
ctx.fillStyle = rb === this.sel ? '#f2c879' : rb.def.color;
ctx.beginPath(); ctx.arc(m[0], m[1], rb === this.sel ? 3.4 : t === 'planet' ? 2.4 : 1.6, 0, TAU); ctx.fill();
this.mmPts.push({ rb, x: m[0] + c, y: m[1] + c });
}
const cp = Cam.pos, cr = Math.hypot(cp[0], cp[1], cp[2]);
if (cr > 1e-9) {
const au = ScaleState.inv(cr), k = au / cr, m = this.mmMap(cp[0] * k, cp[2] * k, R);
const f = V.norm(V.sub(Cam.target, Cam.pos)), ang = Math.atan2(f[2], f[0]);
const mx = clamp(m[0], -R, R), my = clamp(m[1], -R, R);
ctx.fillStyle = 'rgba(143,195,255,0.22)'; ctx.beginPath(); ctx.moveTo(mx, my); ctx.arc(mx, my, 18, ang - 0.45, ang + 0.45); ctx.closePath(); ctx.fill();
ctx.fillStyle = '#e8ecf3'; ctx.beginPath(); ctx.arc(mx, my, 2.6, 0, TAU); ctx.fill();
}
ctx.restore();
},
minimapClick(e) {
const r = e.target.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
let best = null, bd = 12;
(this.mmPts || []).forEach(p => { const d = Math.hypot(p.x - x, p.y - y); if (d < bd) { bd = d; best = p.rb; } });
if (best) this.select(best, { fly: true });
},
toggleGuide(v) {
const g = $('#guide'); g.classList.toggle('open', v); g.setAttribute('aria-hidden', !v);
if (v) { this._guideFrom = document.activeElement; g.querySelector('.guide-card').scrollTop = 0; setTimeout(() => $('#guide-x').focus({ preventScroll: true }), 30); }
else if (this._guideFrom && this._guideFrom.focus) this._guideFrom.focus({ preventScroll: true });
},
toggleHelp(v) { const h = $('#help'); h.classList.toggle('open', v); h.setAttribute('aria-hidden', !v); if (v) $('#help-close').focus(); },
bindKeys() {
document.addEventListener('keydown', e => {
if (Flight.on) return;
const tag = (e.target.tagName || '').toLowerCase();
if (e.key === 'k' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); $('#q').focus(); return; }
if (tag === 'input' || tag === 'textarea' || tag === 'select') return;
if ($('#settings') && $('#settings').classList.contains('open')) return;
const key = Keys.norm(e), act = Keys.action('explore', key);
if (act === 'search') { e.preventDefault(); $('#q').focus(); return; }
if (!this.started) return;
const k = e.key;
if (k === 'Escape') {
if ($('#help').classList.contains('open')) this.toggleHelp(false);
else if (this.tour) this.stopTour();
else if (S.mode === 'cine') this.setMode('explore');
else if (this.section) this.openSection(null);
else if (this.sel) this.select(null);
return;
}
if (e.ctrlKey || e.metaKey || e.altKey) return;
if (k === 'Home') { e.preventDefault(); App.home(); return; }
if (this.tour && (k === ' ' || k === 'ArrowRight' || k === 'ArrowLeft')) {
e.preventDefault();
if (k === ' ') this.pauseTour(!this.tour.paused); else this.stepTour(k === 'ArrowRight' ? 1 : -1);
return;
}
const planets = ['sol', 'mercurio', 'venus', 'tierra', 'marte', 'jupiter', 'saturno', 'urano', 'neptuno', 'pluton'];
if (/^[0-9]$/.test(k)) { this.userActed(); this.select(World.byId[planets[+k]], { fly: true }); return; }
if (k === '?') { this.toggleHelp(true); return; }
const step = 0.08;
if (act && (key === 'space' || key.startsWith('arrow'))) e.preventDefault();
switch (act) {
case 'pause': this.togglePause(); break;
case 'faster': this.speed(1); break;
case 'slower': this.speed(-1); break;
case 'now': Time.goLive(); this.refreshTime(); this.toast(I18N.t('Fecha y hora actuales en tiempo real')); break;
case 'follow': this.camAction('follow'); break;
case 'center': this.camAction('center'); break;
case 'explore': this.camAction('explore'); break;
case 'overview': this.camAction('overview'); break;
case 'sun': this.camAction('sun'); break;
case 'zoomIn': this.camAction('in'); break;
case 'zoomOut': this.camAction('out'); break;
case 'rotL': Cam.cancel(); Cam.az += step; break;
case 'rotR': Cam.cancel(); Cam.az -= step; break;
case 'rotU': Cam.cancel(); Cam.el += step; break;
case 'rotD': Cam.cancel(); Cam.el -= step; break;
case 'names': { const v = !S.layers.namesPlanets; ['namesPlanets', 'namesMoons', 'namesSmall', 'namesCraft'].forEach(x => this.setLayer(x, v)); this.toast(I18N.t(v ? 'Nombres visibles' : 'Nombres ocultos')); break; }
case 'orbits': { const v = !S.layers.orbits; this.setLayer('orbits', v); this.setLayer('moonOrbits', v); this.toast(I18N.t(v ? 'Órbitas visibles' : 'Órbitas ocultas')); break; }
case 'constel': { const v = !S.layers.constLines; this.setLayer('constLines', v); this.setLayer('constNames', v); this.toast(I18N.t(v ? 'Constelaciones visibles' : 'Constelaciones ocultas')); break; }
case 'minimap': this.toggleMinimap(); break;
case 'cine': this.setMode(S.mode === 'cine' ? 'explore' : 'cine'); break;
case 'tours': this.tour ? this.stopTour() : this.toggleTours(true); break;
case 'pilot': Flight.openHangar(); break;
}
});
},
frame(dt) {
this.dateSaveT = (this.dateSaveT || 0) + dt; if (this.dateSaveT > 10) { this.dateSaveT = 0; SettingsApply.saveDate(); }
this.updateLabels();
this.updateLive();
this.tickTour(dt);
if ((World.frame & 3) === 0) { this.tickTime(); this.drawMinimap(); }
},
};
const TOUR_SYSTEMS = [
['tierra', 'la Tierra', ['luna']],
['marte', 'Marte', ['fobos', 'deimos']],
['jupiter', 'Júpiter', ['io', 'europa', 'ganimedes', 'calisto']],
['saturno', 'Saturno', ['mimas', 'encelado', 'tetis', 'dione', 'rea', 'titan', 'japeto']],
['urano', 'Urano', ['miranda', 'ariel', 'umbriel', 'titania', 'oberon']],
['neptuno', 'Neptuno', ['triton']],
['pluton', 'Plutón', ['caronte']],
];
const TOURS = [
{ id: 'planetas', name: 'Planetas', title: 'Planetas del Sistema Solar', desc: 'Explora los ocho planetas principales del Sistema Solar, de Mercurio a Neptuno.',
layers: ['orbits', 'namesPlanets'],
stops: ['mercurio', 'venus', 'tierra', 'marte', 'jupiter', 'saturno', 'urano', 'neptuno'].map(id => ({ id, kind: 'planet', hold: 10 })) },
{ id: 'lunas', name: 'Lunas', title: 'Lunas del Sistema Solar', desc: 'Descubre los satélites naturales más importantes, sistema por sistema.',
layers: ['moons', 'namesMoons', 'moonOrbits'],
stops: TOUR_SYSTEMS.flatMap(([p, art, ms]) => [{ id: p, kind: 'system', hold: 6, group: 'Sistema de ' + art, moons: ms }, ...ms.map(id => ({ id, kind: 'moon', hold: 8, group: 'Sistema de ' + art }))]) },
{ id: 'exploracion', name: 'Exploración humana', title: 'Exploración humana', desc: 'Visita estaciones espaciales, telescopios y sondas que exploran el Sistema Solar.',
layers: ['craft', 'namesCraft'],
stops: [['iss', 'Órbita terrestre baja'], ['tiangong', 'Órbita terrestre baja'], ['hubble', 'Órbita terrestre baja'], ['jwst', 'Puntos de Lagrange'], ['euclid', 'Puntos de Lagrange'], ['soho', 'Puntos de Lagrange'],
['parker', 'Sondas del Sistema Solar'], ['newhorizons', 'Sondas del Sistema Solar'], ['voyager2', 'Sondas del Sistema Solar'], ['voyager1', 'Sondas del Sistema Solar']].map(([id, group]) => ({ id, kind: 'craft', hold: 10, group })) },
{ id: 'menores', name: 'Asteroides y cometas', title: 'Asteroides y cometas', desc: 'Recorre algunos de los objetos menores más relevantes y compara sus trayectorias con las de los planetas.',
layers: ['asteroids', 'namesSmall', 'orbits'],
stops: [['ceres', 'Cinturón de asteroides'], ['vesta', 'Cinturón de asteroides'], ['palas', 'Cinturón de asteroides'], ['higia', 'Cinturón de asteroides'],
['eros', 'Asteroides cercanos a la Tierra'], ['bennu', 'Asteroides cercanos a la Tierra'], ['ryugu', 'Asteroides cercanos a la Tierra'], ['apofis', 'Asteroides cercanos a la Tierra'],
['halley', 'Cometas'], ['halebopp', 'Cometas'], ['c67p', 'Cometas']].map(([id, group]) => ({ id, kind: 'minor', hold: 13, group })) },
{ id: 'clasico', name: 'Gran recorrido', title: 'Gran recorrido del Sistema Solar', desc: 'Un viaje panorámico del Sol al cometa Halley con lo más destacado de cada mundo.',
layers: ['orbits'],
stops: TOUR.map(s => ({ id: s.id, kind: 'classic', hold: 7.5, txt: s.txt })) },
];
const TOUR_ICONS = {
planetas: '<circle cx="12" cy="12" r="4.2"/><ellipse cx="12" cy="12" rx="9.5" ry="3.4" transform="rotate(-20 12 12)"/>',
lunas: '<circle cx="10" cy="13" r="6"/><circle cx="19" cy="6" r="2.2"/>',
exploracion: '<path d="M4 14l4-1 3 3-1 4M9 15l7-7c2-2 4-3 5-3 0 1-1 3-3 5l-7 7"/><circle cx="15.5" cy="8.5" r="1.3"/>',
menores: '<path d="M3 21l7-7"/><circle cx="14" cy="10" r="4.5"/><path d="M5 9c1-2 3-3 5-3"/>',
clasico: '<circle cx="12" cy="12" r="3"/><circle cx="12" cy="12" r="7.5" stroke-dasharray="2 3"/><circle cx="19.5" cy="12" r="1.4"/>',
};
const tourDuration = t => t.stops.reduce((s, x) => s + x.hold + (x.kind === 'minor' ? 6 : 3.5), 0);
Object.assign(UI, {
tour: null, tourRB: null,
buildTours() {
this.renderTourCards();
$('#tours-list').addEventListener('click', e => { const b = e.target.closest('[data-tour]'); if (b) this.launchTour(b.dataset.tour, b.dataset.cine === '1'); });
$('#tours-x').addEventListener('click', () => this.toggleTours(false));
$('#tours').addEventListener('click', e => { if (e.target.id === 'tours') this.toggleTours(false); });
addEventListener('keydown', e => { if (e.key === 'Escape' && $('#tours').classList.contains('open')) { e.stopImmediatePropagation(); this.toggleTours(false); } }, true);
const B = (id, f) => $('#' + id).addEventListener('click', () => { f(); $('#gl').focus({ preventScroll: true }); });
B('tb-prev', () => this.stepTour(-1)); B('tb-next', () => this.stepTour(1)); B('tb-play', () => this.pauseTour(!this.tour.paused));
B('tb-restart', () => this.restartTour()); B('tb-exit', () => this.stopTour()); B('tb-x', () => this.stopTour());
B('tb-card', () => { const rb = this.tourRB; this.stopTour(); if (rb) this.select(rb); });
B('tb-again', () => this.restartTour()); B('tb-other', () => { this.stopTour(); this.toggleTours(true); });
$('#tb-play').innerHTML = ICON.pause;
},
renderTourCards() {
$('#tours-list').innerHTML = TOURS.map(tr => `
<article class="tour-card">
<div class="tour-ico"><svg viewBox="0 0 24 24" aria-hidden="true">${TOUR_ICONS[tr.id]}</svg></div>
<h3>${tr.name}</h3>
<p>${tr.desc}</p>
<p class="tour-meta">${tr.stops.length} paradas · ≈ ${Math.max(1, Math.round(tourDuration(tr) / 60))} min</p>
<div class="tour-act"><button class="cta" data-tour="${tr.id}">Iniciar recorrido</button><button class="cta-ghost tour-cine-btn" data-tour="${tr.id}" data-cine="1"><svg viewBox="0 0 24 24" width="16" height="16" aria-hidden="true"><rect x="3" y="5" width="18" height="14" rx="2" fill="none" stroke="currentColor" stroke-width="1.6"/><path d="M7 5v14M17 5v14M3 9.5h4M3 14.5h4M17 9.5h4M17 14.5h4" fill="none" stroke="currentColor" stroke-width="1.6"/></svg>Iniciar en Modo Cine</button></div>
</article>`).join('');
},
toggleTours(v) {
const g = $('#tours'); g.classList.toggle('open', v); g.setAttribute('aria-hidden', !v);
if (v) { this._toursFrom = document.activeElement; setTimeout(() => { const b = g.querySelector('[data-tour]'); if (b) b.focus({ preventScroll: true }); }, 40); }
else if (this._toursFrom && this._toursFrom.focus) this._toursFrom.focus({ preventScroll: true });
},
launchTour(id, cine) {
this.toggleTours(false);
if (!this.started) App.begin({ quiet: true });
this.startTour(id, { cine });
},
startTour(id, opts) {
opts = opts || {};
const def = TOURS.find(t => t.id === (id || 'clasico'));
if (this.tour) this.stopTour(true);
if (S.mode === 'cine') this.setMode('explore');
(def.layers || []).forEach(l => { if (!S.layers[l]) this.setLayer(l, true); });
this.select(null, { keepCam: true, tour: true }); this.openSection(null);
this.tour = { def, i: -1, t: 0, phase: 'travel', paused: false, done: false };
$('#app').classList.add('touring');
$('#tourbar').classList.add('on'); $('#tourbar').classList.remove('done');
if (opts.cine) {
this.tour.phase = 'intro'; this.enterTourCine();
const T = this.tour; setTimeout(() => { if (this.tour === T) this.gotoStop(0, 1); }, 1400);
} else this.gotoStop(0, 1);
},
gotoStop(i, dir) {
const T = this.tour; if (!T) return;
const n = T.def.stops.length;
while (i >= 0 && i < n && (!World.byId[T.def.stops[i].id] || World.byId[T.def.stops[i].id].hidden)) i += dir || 1;
if (i >= n) return this.finishTour();
if (i < 0) i = T.i >= 0 ? T.i : 0;
T.i = i; T.t = 0; T.phase = 'travel'; T.done = false;
const stop = T.def.stops[i], rb = World.byId[stop.id];
this.tourRB = rb;
let dist;
if (stop.kind === 'system') {
const ext = Math.max(...stop.moons.map(m => (World.byId[m] && World.byId[m].moonDistMean) || 0));
dist = clamp(Math.max(ext * 2.6, Cam.frameDist(rb)), rb.rS * 3, Cam.maxDist() * 0.9);
}
Cam.travel(rb, Object.assign(dist ? { dist } : {}, this.tourCine ? { cine: true } : {}));
Cam.auto = 0; this.hideCineCaption(); T.cam0 = null;
$('#tourbar').classList.remove('done'); $('#tourbar').classList.add('travel');
this.renderTourCard(stop, rb);
this.refreshTourBar(); this.refreshCine();
},
stepTour(d) {
const T = this.tour; if (!T) return;
if (T.done) { if (d < 0) this.gotoStop(T.def.stops.length - 1, -1); return; }
this.gotoStop(T.i + d, d);
},
restartTour() { const T = this.tour; if (!T) return; T.paused = false; this.gotoStop(0, 1); },
pauseTour(v, manual) {
const T = this.tour; if (!T || T.done) return;
T.paused = v; T.manual = !!(v && manual);
if (v) Cam.auto = 0;
this.refreshTourBar(); this.refreshCine();
},
finishTour() {
const T = this.tour; if (!T) return;
this.exitTourCine();          // al terminar, la interfaz vuelve con suavidad y aparece el resumen
T.done = true; T.paused = false; this.tourRB = null;
Cam.travelTo([0, 0, 0], ScaleState.get('overview'), 0.42); Cam.auto = 0;
$('#tourbar').classList.add('done');
$('#tb-done-txt').textContent = 'Has completado «' + T.def.title + '»: ' + T.def.stops.length + ' paradas.';
this.refreshTourBar();
},
stopTour(silent) {
if (!this.tour) return;
this.exitTourCine();
this.tour = null; this.tourRB = null; Cam.auto = S.mode === 'cine' ? 0.035 : 0;
$('#tourbar').classList.remove('on', 'done', 'travel'); $('#app').classList.remove('touring');
if (!silent && this.section === 'explorar') this.openSection('explorar');
},
camTouched() { if (this.tour && !this.tour.done && !this.tour.paused) this.pauseTour(true, true); else if (!this.tour) this.userActed(); },
tickTour(dt) {
const T = this.tour; if (!T || T.done || T.i < 0 || T.phase === 'intro') return;
const stop = T.def.stops[T.i], rb = World.byId[stop.id];
const hold = stop.hold * (this.tourCine ? 1.35 : 1);       // en Modo Cine se contempla cada objeto más tiempo
const prog = k => { $('#tb-prog').style.transform = `scaleX(${k})`; const c = $('#cc-prog'); if (c) c.style.transform = `scaleX(${k})`; };
if (T.paused) { prog(clamp(T.t / hold, 0, 1)); return; }
if (T.phase === 'travel') { if (!Cam.fly) { T.phase = 'hold'; T.t = 0; $('#tourbar').classList.remove('travel'); if (this.tourCine) this.cineCaption(stop, rb); } return; }
T.t += dt;
Cam.auto = T.phase === 'orbit' ? 0.02 : 0.05;
if (this.tourCine && T.phase === 'hold' && !Cam.fly && !UI.reducedMotion) {
if (!T.cam0) T.cam0 = { el: Cam.el, d: Cam.dDist };
const k = T.t;
Cam.auto = 0.032 + 0.018 * Math.sin(k * 0.21);
Cam.el = clamp(T.cam0.el + 0.15 * Math.sin(k * 0.17), -1.2, 1.2);
Cam.dDist = T.cam0.d * (1 - 0.16 * smoothstep(0, hold, k));
}
if (stop.kind === 'minor' && T.phase === 'hold' && T.t > 5 && rb.orbitBuf && !this.tourCine) {
const b = rb.orbitBuf, lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
for (let k = 0; k < b.length; k += 3) for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], b[k + a]); hi[a] = Math.max(hi[a], b[k + a]); }
const c = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2], r = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2;
Cam.travelTo(c, clamp(r * 2.3, rb.rS * 10, Cam.maxDist() * 0.95), 0.95);
T.phase = 'orbit';
$('#tb-phase').textContent = 'Trayectoria completa alrededor del Sol';
}
if (T.phase === 'orbit' && Cam.fly) return;
prog(clamp(T.t / hold, 0, 1));
if (T.t >= hold) this.gotoStop(T.i + 1, 1);
T.live = (T.live || 0) + dt; if (T.live > 0.5) { T.live = 0; this.updateTourLive(); }
},
renderTourCard(stop, rb) {
const T = this.tour, d = rb.def, n = T.def.stops.length;
$('#tb-tour').textContent = 'Recorrido: ' + T.def.title + (stop.group ? ' · ' + stop.group : '');
$('#tb-step').textContent = 'Parada ' + (T.i + 1) + ' de ' + n + ' — ' + (stop.kind === 'system' ? stop.group : d.name);
$('#tb-phase').textContent = '';
let kicker = TYPE_LABEL[d.type] || '', title = d.name, text = d.info.desc, note = '';
if (d.type === 'moon') kicker = 'Luna de ' + rb.parent.def.name;
if (rb.isCraft) kicker = d.sub;
if (stop.kind === 'classic') text = stop.txt;
if (stop.kind === 'system') {
const names = stop.moons.map(m => BODY[m].name), cnt = d.info.moons;
title = stop.group.charAt(0).toUpperCase() + stop.group.slice(1); kicker = TYPE_LABEL[d.type];
const art = TOUR_SYSTEMS.find(s => s[0] === stop.id)[1];
const who = I18N.lang === 'es' ? art.charAt(0).toUpperCase() + art.slice(1) : d.name;
const list = names.length > 1 ? names.slice(0, -1).join(', ') + ' ' + t('y') + ' ' + names[names.length - 1] : names[0];
text = cnt === 1 ? t('{0} tiene una sola luna: {1}.', who, names[0])
: t('{0} tiene {1} lunas conocidas; en este recorrido visitamos {2}: {3}.', who, cnt, names.length === 1 ? t('una') : names.length, list);
}
if (stop.kind === 'minor') note = this.minorNote(rb);
$('#tb-kicker').textContent = kicker; $('#tb-title').textContent = title; $('#tb-text').textContent = text;
$('#tb-note').textContent = note; $('#tb-note').hidden = !note;
$('#tb-data').innerHTML = this.tourFields(stop, rb).map(f => `<div><dt>${esc(f.label)}</dt><dd>${esc(f.val)}</dd></div>`).join('');
this.updateTourLive(true);
},
tourFields(stop, rb) {
if (rb.isCraft) {
const c = rb.def.craft;
return [{ label: 'Tipo de misión', val: rb.def.sub }, { label: 'Agencia', val: c.agency }, { label: 'Lanzamiento', val: c.launchTxt },
{ label: 'Estado', val: c.status }, { label: 'Ubicación', val: c.orbitTxt }];
}
const F = Info.fields(rb), all = [...F.phys, ...F.orb, ...F.comp];
const pick = labels => labels.map(l => all.find(f => typeof l === 'string' ? f.label === l : l.test(f.label))).filter(f => f && f.val != null);
if (stop.kind === 'system') return pick(['Número de lunas', 'Diámetro medio']);
if (rb.def.type === 'moon') return pick(['Diámetro medio', /^Periodo orbital/, /^Distancia media a/, 'Gravedad superficial', 'Temperatura aproximada']);
if (stop.kind === 'minor') return pick(['Diámetro medio', 'Duración de un año (periodo orbital)', 'Excentricidad', 'Perihelio y afelio', 'Inclinación orbital']);
return pick(['Diámetro medio', 'Gravedad superficial', 'Duración de un día', 'Duración de un año (periodo orbital)', 'Número de lunas', 'Temperatura aproximada']);
},
minorNote(rb) {
const o = rb.def.orbit; if (!o || !(o.t === 'jpl' || o.t === 'kep')) return '';
const k = Astro.elementsAt(o, World.jd), em = Astro.elementsAt(World.byId.mercurio.def.orbit, World.jd).e;
const q = k.a * (1 - k.e), Q = k.a * (1 + k.e);
let s = k.e > 0.5 ? `Su órbita es muy alargada (excentricidad ${fmt(k.e, 3)}); entre los planetas, la más excéntrica es la de Mercurio (${fmt(em, 3)}).`
: k.e > em ? `Su órbita es más excéntrica que la de cualquier planeta (Mercurio: ${fmt(em, 3)}).`
: k.e > 0.1 ? `Su órbita es algo alargada (excentricidad ${fmt(k.e, 3)}), comparable a la de Mercurio (${fmt(em, 3)}), el planeta de órbita más excéntrica.`
: `Su órbita es casi circular (excentricidad ${fmt(k.e, 3)}), como la de la mayoría de los planetas.`;
if (k.i > 90) s += ' Además, recorre su órbita en sentido retrógrado, al revés que los planetas.';
else if (k.i > 10) s += ` Está inclinada ${fmt(k.i, 1)}° respecto al plano en que giran los planetas.`;
if (q < 1.017 && Q > 0.983 && k.i < 30) s += ' Su órbita cruza la de la Tierra.';
else if (q < 0.983) s += ' En su perihelio se acerca al Sol más que la Tierra.';
else if (q < 1.3) s += ' Se acerca a la órbita terrestre: es un objeto cercano a la Tierra.';
return s;
},
updateTourLive(force) {
const rb = this.tourRB, el = $('#tb-live'); if (!rb || !el) return;
const rows = Info.live(rb).filter(r => /^Distancia actual|^Altitud|^Tiempo de luz/.test(r.label)).slice(0, 2);
el.innerHTML = rows.map(r => `<span>${esc(r.label)}: <b>${r.txt ? esc(r.txt) : esc(fmtKm(r.km))}</b>${r.au != null && r.au >= 0.01 ? ` <small>(${fmtAU(r.au)})</small>` : ''}</span>`).join('');
},
refreshTourBar() {
const T = this.tour; if (!T) return;
const p = $('#tb-play');
p.innerHTML = T.paused ? ICON.play : ICON.pause;
p.setAttribute('aria-label', T.paused ? 'Reanudar recorrido' : 'Pausar recorrido');
p.dataset.tip = (T.paused ? 'Reanudar' : 'Pausar') + ' (espacio)';
$('#tb-paused').textContent = T.paused ? (T.manual ? 'En pausa: mueve la cámara libremente y pulsa Reanudar para continuar' : 'Recorrido en pausa') : '';
$('#tb-prev').disabled = !T.done && T.i <= 0;
$('#tb-dots').style.setProperty('--p', T.done ? 1 : (T.i + 1) / T.def.stops.length);
},
});
Object.assign(UI, {
initAbout() {
$('#about-btn').addEventListener('click', () => this.toggleAbout(true));
$('#about-x').addEventListener('click', () => this.toggleAbout(false));
$('#about-ok').addEventListener('click', () => this.toggleAbout(false));
$('#about').addEventListener('click', e => { if (e.target.id === 'about') this.toggleAbout(false); });
addEventListener('keydown', e => { if (e.key === 'Escape' && $('#about').classList.contains('open')) { e.stopImmediatePropagation(); this.toggleAbout(false); } }, true);
},
toggleAbout(v) {
const g = $('#about'); g.classList.toggle('open', v); g.setAttribute('aria-hidden', !v);
if (v) { this._aboutFrom = document.activeElement; g.querySelector('.about-card').scrollTop = 0; setTimeout(() => $('#about-x').focus({ preventScroll: true }), 30); }
else { if (this._aboutFrom && this._aboutFrom.focus) this._aboutFrom.focus({ preventScroll: true }); }
},
});
Object.assign(UI, {
creditAssets() {
const site = url => { try { const h = new URL(url).hostname.replace(/^www\./, ''); return h === 'sketchfab.com' ? 'Sketchfab' : h; } catch (e) { return null; } };
const fromCredit = (use, c) => ({ use, name: c.title, author: c.author, source: c.url ? site(c.url) : null, license: c.license, licenseUrl: c.licenseUrl, url: c.url, note: c.note });
const groups = [
['Telescopios, estaciones y sondas', [['Telescopio espacial James Webb', 'jwst'], ['Telescopio espacial Hubble', 'hubble'], ['Voyager 1 y Voyager 2', 'voyager'], ['Estación Espacial Internacional', 'iss'], ['Estación espacial Tiangong', 'tiangong'], ['Sonda solar Parker', 'parker'], ['Sonda New Horizons', 'newhorizons']]
.filter(([, k]) => MODEL_CREDITS[k]).map(([u, k]) => fromCredit(u, MODEL_CREDITS[k]))],
['El Sol, planetas y lunas', Object.keys(PLANET_TEX).map(id => fromCredit(BODY[id].name, BODY[id].credit))],
['Naves del modo de vuelo', SHIPS.filter(sh => MODEL_CREDITS[sh.id]).map(sh => fromCredit('Nave ' + sh.name, MODEL_CREDITS[sh.id]))],
['Recursos utilizados en versiones anteriores', PREV_CREDITS.map(p => fromCredit(BODY[p.id].name + ' (versión anterior)', p.credit))],
['Música', [{ use: 'Tema musical principal de SOLARIS', name: MUSIC_INFO.title, author: MUSIC_INFO.artist, source: 'Archivo aportado por el creador del proyecto', license: null, url: null }]],
['Referencias culturales', SHIPS.filter(sh => sh.brand && sh.brand.credit).map(sh => ({ use: 'Ficha de la nave ' + sh.name, name: sh.brand.credit.name, author: sh.brand.credit.author, source: 'Imagen aportada por el creador del proyecto',
license: sh.id === 'ranger' ? 'Marca y obra protegidas; se muestran como referencia cultural (ver «Referencias a Interstellar»)' : 'Información de atribución pendiente; marca de sus titulares, mostrada como referencia cultural (ver «Referencias a ' + sh.brand.alt + '»)', url: null }))],
];
return groups.filter(g => g[1].length);
},
buildCredits() {
document.querySelectorAll('#credits [data-ship-ref]').forEach(el => { el.hidden = !SHIPS.some(s => s.id === el.dataset.shipRef); });
const na = '<span class="cr-na">No indicado</span>';
$('#cr-assets').innerHTML = this.creditAssets().map(([title, items]) => `
<h4>${esc(title)}</h4>
<div class="cr-asset-grid">${items.map(a => `
<dl class="cr-asset">
<div class="cr-use">${esc(a.use)}</div>
<div><dt>Nombre del recurso</dt><dd>${esc(a.name || '—')}</dd></div>
<div><dt>Autor</dt><dd>${a.author ? esc(a.author) : na}</dd></div>
<div><dt>Fuente</dt><dd>${a.source ? esc(a.source) : 'Archivo aportado por el creador del proyecto'}</dd></div>
<div><dt>Licencia</dt><dd>${a.license ? (a.licenseUrl ? `<a href="${a.licenseUrl}" target="_blank" rel="noopener">${esc(a.license)}</a>` : esc(a.license)) : na}</dd></div>
<div><dt>Enlace original</dt><dd>${a.url ? `<a href="${a.url}" target="_blank" rel="noopener">${esc(a.url.replace(/^https?:\/\//, ''))}</a>` : '<span class="cr-na">No disponible</span>'}</dd></div>
${a.note ? `<div class="cr-adapt">${esc(a.note.charAt(0).toUpperCase() + a.note.slice(1))}.</div>` : ''}
</dl>`).join('')}</div>`).join('');
const back = () => this.closeCredits();
$('#cr-back').addEventListener('click', back); $('#cr-back2').addEventListener('click', back);
$('#credits-btn').addEventListener('click', () => this.openCredits());
$('#about-credits').addEventListener('click', () => { this.toggleAbout(false); this.openCredits(); });
document.querySelectorAll('.cr-index a').forEach(a => a.addEventListener('click', e => {
e.preventDefault();
const t = document.querySelector(a.getAttribute('href')), sc = $('#cr-scroll');
if (t) sc.scrollTo({ top: t.getBoundingClientRect().top - sc.getBoundingClientRect().top + sc.scrollTop - 84, behavior: this.reducedMotion ? 'auto' : 'smooth' });
}));
addEventListener('keydown', e => { if (e.key === 'Escape' && $('#credits').classList.contains('open')) { e.stopImmediatePropagation(); this.closeCredits(); } }, true);
addEventListener('popstate', () => {
if (location.hash === '#creditos') this.openCredits(false);
else if ($('#credits').classList.contains('open')) this.closeCredits(true);
});
if (location.hash === '#creditos') setTimeout(() => this.openCredits(false), 50);
},
openCredits(push) {
const el = $('#credits');
if (el.classList.contains('open')) return;
this._crFrom = document.activeElement; this._crFromIntro = !this.started;
if (push !== false && location.hash !== '#creditos') { try { history.pushState({ solaris: 'creditos' }, '', '#creditos'); } catch (e) {} }
$('#cr-scroll').scrollTop = 0;
const els = [...document.querySelectorAll('#credits .cr-in')];
els.forEach(x => x.classList.remove('vis'));
els.forEach((x, i) => setTimeout(() => x.classList.add('vis'), this.reducedMotion ? 0 : 120 + Math.min(i, 5) * 90));
el.classList.add('open'); el.setAttribute('aria-hidden', 'false'); $('#app').classList.add('in-credits');
setTimeout(() => $('#cr-back').focus({ preventScroll: true }), 60);
},
closeCredits(fromHistory) {
const el = $('#credits');
if (!el.classList.contains('open')) return;
if (!fromHistory && history.state && history.state.solaris === 'creditos') { history.back(); return; }
el.classList.remove('open'); el.setAttribute('aria-hidden', 'true'); $('#app').classList.remove('in-credits');
if (location.hash === '#creditos') { try { history.replaceState(null, '', location.pathname + location.search); } catch (e) {} }
if (this.started && this._crFromIntro) App.home();
else if (this._crFrom && this._crFrom.focus) this._crFrom.focus({ preventScroll: true });
},
});
const SFX = {
ctx: null, out: null, bus: {}, noise: null,
get prefs() { return Settings.state.audio; },            // settings.audio es la única fuente
last: {}, voices: {}, cruise: null, hoverEl: null, hoverT: 0,
ensure() {
if (this.ctx) return this.ctx;
const AC = window.AudioContext || window.webkitAudioContext; if (!AC) return null;
try { this.ctx = new AC({ latencyHint: 'interactive' }); } catch (e) { this.ctx = new AC(); }
const ctx = this.ctx;
const comp = ctx.createDynamicsCompressor();
comp.threshold.value = -16; comp.knee.value = 12; comp.ratio.value = 3; comp.attack.value = 0.004; comp.release.value = 0.25;
this.out = ctx.createGain(); this.out.connect(comp); comp.connect(ctx.destination);
['music', 'ui', 'navigation', 'spacecraft', 'alerts'].forEach(k => { const g = ctx.createGain(); g.connect(this.out); this.bus[k] = g; });
const len = ctx.sampleRate, buf = ctx.createBuffer(1, len, ctx.sampleRate), d = buf.getChannelData(0);
for (let i = 0; i < len; i++) d[i] = Math.random() * 2 - 1;
this.noise = buf;
this.apply(0);
return ctx;
},
apply(tc) {
if (!this.ctx) return;
const t = this.ctx.currentTime, k = tc == null ? 0.08 : tc, P = this.prefs;
const set = (g, v) => { g.gain.cancelScheduledValues(t); if (k) g.gain.setTargetAtTime(v, t, k); else g.gain.value = v; };
set(this.out, P.muted ? 0 : P.master);
set(this.bus.music, 1);
const cine = UI.tourCine ? 0.3 : 1;     // Modo Cine: la interfaz no rompe la inmersión
set(this.bus.ui, P.ui * 2 * cine); set(this.bus.navigation, P.nav * 1.9);
set(this.bus.spacecraft, P.ship); set(this.bus.alerts, P.ship * 0.9);
},
set(key, v) { Settings.set('audio.' + key, v); },
unlock() { const c = this.ensure(); if (c && c.state === 'suspended') c.resume(); },
ready() { return this.ctx && this.ctx.state === 'running' && !this.prefs.muted; },
tone(bus, o) {
const c = this.ctx, t = c.currentTime + (o.delay || 0), dur = o.dur, nodes = [];
const osc = c.createOscillator(), g = c.createGain();
osc.type = o.type || 'sine'; osc.frequency.setValueAtTime(o.f, t);
if (o.f2) osc.frequency.exponentialRampToValueAtTime(o.f2, t + dur * (o.glide || 0.85));
if (o.detune) osc.detune.value = o.detune;
g.gain.setValueAtTime(0.0001, t); g.gain.linearRampToValueAtTime(o.gain, t + (o.a || 0.004));
g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
let last = g;
osc.connect(g);
if (o.lp) { const f = c.createBiquadFilter(); f.type = 'lowpass'; f.frequency.value = o.lp; g.connect(f); last = f; nodes.push(f); }
last.connect(this.bus[bus]); osc.start(t); osc.stop(t + dur + 0.03);
nodes.push(osc, g);
return { nodes, g, end: t + dur };
},
hiss(bus, o) {
const c = this.ctx, t = c.currentTime + (o.delay || 0), dur = o.dur;
const s = c.createBufferSource(), f = c.createBiquadFilter(), g = c.createGain();
s.buffer = this.noise; s.loop = true; s.playbackRate.value = 0.7 + Math.random() * 0.6;
f.type = 'bandpass'; f.Q.value = o.q || 1; f.frequency.setValueAtTime(o.f0, t);
if (o.f1) f.frequency.exponentialRampToValueAtTime(o.f1, t + dur * (o.sweep || 0.9));
const pk = o.peak == null ? 0.35 : o.peak;
g.gain.setValueAtTime(0.0001, t); g.gain.linearRampToValueAtTime(o.gain, t + dur * pk);
g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
s.connect(f); f.connect(g); g.connect(this.bus[bus]); s.start(t); s.stop(t + dur + 0.03);
return { nodes: [s, f, g], g, end: t + dur };
},
SOUNDS: {
hover: { bus: 'ui', gap: 0.07, max: 1, fn: S => [S.tone('ui', { f: 740, dur: 0.035, gain: 0.016, lp: 2200 }), S.hiss('ui', { f0: 2600, dur: 0.012, gain: 0.006, q: 2, peak: 0.2 })] },
click: { bus: 'ui', gap: 0.05, max: 2, fn: S => [S.tone('ui', { f: 587, dur: 0.1, gain: 0.07, lp: 3000 }), S.tone('ui', { f: 880, dur: 0.09, gain: 0.04, delay: 0.022, lp: 3000 }), S.hiss('ui', { f0: 3200, dur: 0.014, gain: 0.018, q: 1.5, peak: 0.15 })] },
minor: { bus: 'ui', gap: 0.05, max: 2, fn: S => [S.tone('ui', { f: 523, dur: 0.055, gain: 0.045, lp: 2400 }), S.hiss('ui', { f0: 2400, dur: 0.009, gain: 0.012, q: 1.6, peak: 0.15 })] },
back: { bus: 'ui', gap: 0.06, max: 2, fn: S => [S.tone('ui', { f: 622, f2: 415, dur: 0.12, gain: 0.06, lp: 2600 }), S.tone('ui', { f: 466, f2: 311, dur: 0.12, gain: 0.03, delay: 0.035, lp: 2000 })] },
on: { bus: 'ui', gap: 0.06, max: 2, fn: S => [S.tone('ui', { f: 440, dur: 0.07, gain: 0.05, lp: 2600 }), S.tone('ui', { f: 659, dur: 0.1, gain: 0.055, delay: 0.06, lp: 2600 })] },
off: { bus: 'ui', gap: 0.06, max: 2, fn: S => [S.tone('ui', { f: 659, dur: 0.07, gain: 0.05, lp: 2600 }), S.tone('ui', { f: 440, dur: 0.11, gain: 0.045, delay: 0.06, lp: 2000 })] },
select: { bus: 'navigation', gap: 0.12, max: 1, fn: S => [
S.tone('navigation', { type: 'triangle', f: 330, f2: 494, dur: 0.16, gain: 0.05, lp: 1600, glide: 0.9 }),
S.tone('navigation', { f: 494, dur: 0.55, gain: 0.04, delay: 0.12, a: 0.02, detune: -4 }),
S.tone('navigation', { f: 740, dur: 0.6, gain: 0.03, delay: 0.15, a: 0.03, detune: 5 }),
S.hiss('navigation', { f0: 1800, f1: 3400, dur: 0.18, gain: 0.012, q: 3, delay: 0.1 })] },
arrive: { bus: 'navigation', gap: 0.3, max: 1, fn: (S, d) => [
S.hiss('navigation', { f0: 1300 - d * 400, f1: 260, dur: 0.45 + d * 0.25, gain: 0.028 + d * 0.012, q: 0.9, peak: 0.08, sweep: 0.8 }),
S.tone('navigation', { f: 587, dur: 0.7, gain: 0.04, delay: 0.16 + d * 0.06, a: 0.015 }),
S.tone('navigation', { f: 880, dur: 0.85, gain: 0.03, delay: 0.25 + d * 0.06, a: 0.02, detune: 3 })] },
},
play(name, arg) {
const def = this.SOUNDS[name]; if (!def || !this.ready()) return;
const now = this.ctx.currentTime;
if (now - (this.last[name] || -9) < def.gap) return;            // sin repeticiones acumuladas
this.last[name] = now;
const list = (this.voices[name] = (this.voices[name] || []).filter(v => v.end > now));
while (list.length >= def.max) this.kill(list.shift(), 0.02);      // reutiliza la voz más antigua
const parts = def.fn(this, arg || 0);
list.push({ parts, end: Math.max(...parts.map(p => p.end)) });
},
kill(v, tc) {
const t = this.ctx.currentTime;
(v.parts || [v]).forEach(p => { try { p.g.gain.cancelScheduledValues(t); p.g.gain.setTargetAtTime(0.0001, t, tc || 0.02); p.nodes.forEach(n => n.stop && n.stop(t + 0.12)); } catch (e) {} });
},
navStart(fly, depth) {
if (!this.ready() || !fly) return;
const d = clamp(depth, 0, 1), dur = fly.dur;
this.navStop(0.08);
const parts = [
this.tone('navigation', { f: 196, f2: 262, dur: 0.22, gain: 0.03 + d * 0.01, lp: 1200, delay: 0.04 }),
this.hiss('navigation', { f0: 240, f1: 900 + d * 900, dur: Math.min(dur * 0.55, 0.5 + d * 0.6), gain: 0.03 + d * 0.025, q: 0.8, peak: 0.45, sweep: 1 }),
];
if (d > 0.35) parts.push(this.tone('navigation', { f: 62, f2: 48, dur: 0.9 + d * 0.6, gain: 0.05 * d, a: 0.25 }));
let cruise = null;
if (dur > 1.6) {
const c = this.ctx, t = c.currentTime;
const o1 = c.createOscillator(), o2 = c.createOscillator(), n = c.createBufferSource(), nf = c.createBiquadFilter(), lp = c.createBiquadFilter(), g = c.createGain();
o1.type = o2.type = 'sawtooth'; o1.frequency.value = 55 - d * 7; o2.frequency.value = (55 - d * 7) * 1.006;
n.buffer = this.noise; n.loop = true; nf.type = 'bandpass'; nf.frequency.value = 420; nf.Q.value = 0.7;
lp.type = 'lowpass'; lp.frequency.value = 140; lp.Q.value = 0.4; g.gain.value = 0.0001;
o1.connect(lp); o2.connect(lp); n.connect(nf); nf.connect(lp); lp.connect(g); g.connect(this.bus.navigation);
o1.start(t); o2.start(t); n.start(t);
cruise = { nodes: [o1, o2, n], lp, nf, g, d };
}
this.cruise = { cruise, parts, d, arrived: false };
},
navUpdate(k) {
const C = this.cruise; if (!C || !C.cruise || !this.ctx) return;
const s = Math.sin(Math.PI * clamp(k, 0, 1)), t = this.ctx.currentTime, X = C.cruise;
X.g.gain.setTargetAtTime(0.0001 + s * (0.022 + X.d * 0.02), t, 0.12);
X.lp.frequency.setTargetAtTime(120 + s * (380 + X.d * 300), t, 0.12);
X.nf.frequency.setTargetAtTime(300 + s * 700, t, 0.15);
},
navArrive() {
const C = this.cruise; if (!C || C.arrived) return;
C.arrived = true; this.navStop(0.25, true);
this.play('arrive', C.d);
},
navStop(tc, keepParts) {
const C = this.cruise; if (!C || !this.ctx) return;
const t = this.ctx.currentTime;
if (C.cruise) { C.cruise.g.gain.cancelScheduledValues(t); C.cruise.g.gain.setTargetAtTime(0.0001, t, tc || 0.08); C.cruise.nodes.forEach(n => { try { n.stop(t + (tc || 0.08) * 6); } catch (e) {} }); C.cruise = null; }
if (!keepParts) C.parts.forEach(p => this.kill(p, 0.05));
if (!keepParts) this.cruise = null;
},
navCancel() { if (this.cruise) { this.navStop(0.1); } },
classify(el) {
if (el.dataset.sfx) return el.dataset.sfx;
if (el.matches('.guide-x, .tb-x, #tb-exit, #cr-back, #cr-back2, #hg-back, #hg-home, #help-close, #about-ok, #guide-ok, #info-close, [aria-label^="Cerrar"], [aria-label^="Salir"]') || /^(←\s*)?(Volver|Cerrar|Salir)/.test((el.textContent || '').trim())) return 'back';
if (el.hasAttribute('aria-pressed')) return el.getAttribute('aria-pressed') === 'true' ? 'off' : 'on';
if (el.matches('.cta, .cta-ghost, [data-tour], .hg-ship, #fl-launch, .pill, .intro-about, #credits-btn, .tour-card button, #tb-again, #tb-other')) return 'click';
return 'minor';
},
init() {
const unlock = () => this.unlock();
['pointerdown', 'keydown', 'touchstart'].forEach(t => addEventListener(t, unlock, { capture: true, passive: true }));
document.addEventListener('click', e => {
const el = e.target.closest('button, a[href], [role="button"], [role="option"]');
if (!el || el.disabled || el.closest('#gl')) return;
if (el.matches('[data-go]') || el.closest('[data-go]')) return;   // la selección de objetos tiene su propio sonido
this.unlock(); this.play(this.classify(el));
}, true);
document.addEventListener('change', e => {
const el = e.target;
if (el.matches('input[type="checkbox"]')) this.play(el.checked ? 'on' : 'off');
else if (el.matches('select, input[type="date"], input[type="datetime-local"]')) this.play('minor');
}, true);
const HOVER = '.cta, .cta-ghost, .intro-mini, .intro-about, .intro-music, .tour-card, .hg-ship, #modes button, #rail button, .fly-btn, .pill, .row, .tb-btn';
document.addEventListener('pointerover', e => {
if (e.pointerType !== 'mouse' || !this.prefs.hover) return;
const el = e.target.closest(HOVER); if (!el || el === this.hoverEl || (e.relatedTarget && el.contains(e.relatedTarget))) return;
this.hoverEl = el; this.hover();
}, true);
document.addEventListener('pointerout', e => { const el = e.target.closest(HOVER); if (el && el === this.hoverEl && !(e.relatedTarget && el.contains(e.relatedTarget))) this.hoverEl = null; }, true);
},
hover() {
if (!this.prefs.hover || UI.tourCine) return;
const now = performance.now(); if (now - this.hoverT < 90) return;
this.hoverT = now; this.play('hover');
},
};
const MIX_ROWS = [['master', 'Volumen general'], ['music', 'Música'], ['ui', 'Efectos de interfaz'], ['nav', 'Navegación'], ['ship', 'Naves y motores']];
const mixerHTML = () => `
<div class="mix">
${MIX_ROWS.map(([k, l]) => `<label class="mix-row"><span>${l}</span><input type="range" min="0" max="100" step="1" data-mix="${k}" aria-label="${l}"><output data-mix-val="${k}"></output></label>`).join('')}
<label class="mix-sw"><input type="checkbox" data-mix-sw="hover"><span class="sw" aria-hidden="true"></span><span>Sonidos al pasar el cursor</span></label>
<label class="mix-sw"><input type="checkbox" data-mix-sw="muted"><span class="sw" aria-hidden="true"></span><span>Silenciar todos los sonidos</span></label>
</div>`;
Object.assign(UI, {
initAudio() {
SFX.init();
document.addEventListener('input', e => {
const k = e.target.dataset && e.target.dataset.mix; if (!k) return;
const v = +e.target.value / 100;
Settings.set('audio.' + k, v);
});
document.addEventListener('change', e => {
const k = e.target.dataset && e.target.dataset.mixSw; if (!k) return;
SFX.set(k, e.target.checked);
});
},
refreshAudio() {
const val = k => Settings.state.audio[k];
document.querySelectorAll('[data-mix]').forEach(i => { if (document.activeElement !== i) i.value = Math.round(val(i.dataset.mix) * 100); });
document.querySelectorAll('[data-mix-val]').forEach(o => { o.textContent = Math.round(val(o.dataset.mixVal) * 100) + ' %'; });
document.querySelectorAll('[data-mix-sw]').forEach(i => { i.checked = !!SFX.prefs[i.dataset.mixSw]; });
document.querySelectorAll('.mix').forEach(m => m.classList.toggle('all-muted', SFX.prefs.muted));
document.querySelectorAll('[data-music]').forEach(b => b.classList.toggle('all-muted', SFX.prefs.muted));
},
});
const MUSIC_INFO = { title: 'Metamorphosis', artist: 'Laura Platt', role: 'Tema musical principal' };
const Music = {
el: null, url: null, ctx: null, gain: null, ready: false, playing: false, duck: 1, waiting: false, fadeTimer: null,
get prefs() { return Settings.state.audio; },
ensureSource() {
if (this.el) return true;
if (ASSET_MANIFEST.music) {           // versión web: la pista es un archivo que se descarga al empezar a sonar
const el = this.el = new Audio(), M = ASSET_MANIFEST.music;
el.src = el.canPlayType('audio/webm; codecs=opus') ? M.webm : M.mp3;
el.loop = true; el.preload = 'auto'; el.crossOrigin = 'anonymous';
return true;
}
const node = document.getElementById('solaris-music'); if (!node) return false;
try {
const txt = node.textContent.trim(), b64 = txt.slice(txt.indexOf(',') + 1), bin = atob(b64), u8 = new Uint8Array(bin.length);
for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
this.url = URL.createObjectURL(new Blob([u8], { type: 'audio/webm' }));
node.textContent = '';                      // libera la copia de texto
const el = this.el = new Audio();
el.src = this.url; el.loop = true; el.preload = 'auto'; el.crossOrigin = 'anonymous';
el.addEventListener('ended', () => { if (!el.loop) el.play(); });
if (!el.canPlayType('audio/webm; codecs=opus')) { this.unsupported = true; UI.refreshMusic && UI.refreshMusic(); }
return true;
} catch (e) { console.warn('No se pudo preparar la música', e); return false; }
},
ensureGraph() {
if (this.ctx || this.noGraph) return;
const ctx = SFX.ensure();
if (!ctx) { this.noGraph = true; return; }
try {
this.ctx = ctx;
const src = ctx.createMediaElementSource(this.el);
this.gain = ctx.createGain(); this.gain.gain.value = 0;
src.connect(this.gain); this.gain.connect(SFX.bus.music);   // bus «music» del sistema global
} catch (e) { this.noGraph = true; this.ctx = null; }
},
level() { return this.prefs.music * this.duck; },
applyLevel() { if (this.playing) this.fadeTo(this.level(), 0.15); },
curGain(now) {
const r = this.ramp; if (!r) return 0;
if (now >= r.t1) return r.v1;
return r.v0 + (r.v1 - r.v0) * Math.max(0, (now - r.t0) / (r.t1 - r.t0));
},
fadeTo(v, sec, done) {
clearTimeout(this.fadeTimer);
if (this.gain) {
const g = this.gain.gain, now = this.ctx.currentTime, cur = this.curGain(now), end = now + Math.max(0.02, sec);
g.cancelScheduledValues(now); g.setValueAtTime(cur, now); g.linearRampToValueAtTime(Math.max(0, v), end);
this.ramp = { v0: cur, v1: Math.max(0, v), t0: now, t1: end };
} else if (this.el) {
const from = this.el.volume, t0 = performance.now(), dur = Math.max(20, sec * 1000);
const step = () => { const k = Math.min(1, (performance.now() - t0) / dur); this.el.volume = clamp(from + (v - from) * k, 0, 1); if (k < 1) this.fadeTimer = setTimeout(step, 30); };
step();
}
if (done) this.fadeTimer = setTimeout(done, sec * 1000 + 40);
},
async play(fade) {
if (!this.prefs.musicOn || !this.ensureSource() || this.unsupported) return;
this.ensureGraph();
try {
if (this.ctx && this.ctx.state !== 'running') await Promise.race([this.ctx.resume(), new Promise(r => setTimeout(r, 400))]);
if (this.ctx && this.ctx.state !== 'running') throw new Error('audio bloqueado hasta la primera interacción');
if (!this.gain) this.el.volume = 0;
await this.el.play();
this.playing = true; this.waiting = false;
this.fadeTo(this.level(), fade == null ? 3.5 : fade);
} catch (e) { this.playing = false; this.waitForGesture(); }
UI.refreshMusic();
},
pause(fade) {
if (!this.el || !this.playing) { this.playing = false; UI.refreshMusic(); return; }
this.playing = false; UI.refreshMusic();
this.fadeTo(0, fade == null ? 1.4 : fade, () => { if (!this.playing) this.el.pause(); });
},
toggle() {
const on = !(this.prefs.musicOn && (this.playing || this.waiting));
if (!on) this.waiting = false;
if (on === this.prefs.musicOn) { if (on) this.play(2.2); else this.pause(); }
Settings.set('audio.musicOn', on);
UI.refreshMusic();
},
setVolume(v) { Settings.set('audio.music', clamp(v, 0, 1)); },
setDuck(f) {
if (Math.abs(f - this.duck) < 1e-3) return;
this.duck = f; if (this.playing) this.fadeTo(this.level(), 1.0);
},
waitForGesture() {
if (this.waiting || !this.prefs.musicOn) return;
this.waiting = true;
const go = () => { ['pointerdown', 'keydown', 'touchend'].forEach(t => removeEventListener(t, go, true)); if (this.waiting && this.prefs.musicOn) { this.waiting = false; this.play(3.5); } };
['pointerdown', 'keydown', 'touchend'].forEach(t => addEventListener(t, go, true));
},
init() {
document.addEventListener('visibilitychange', () => {
if (document.hidden) { if (this.playing) { this.resumeOnShow = true; this.pause(0.5); } }
else if (this.resumeOnShow) { this.resumeOnShow = false; this.play(1.8); }
});
setTimeout(() => { if (this.prefs.musicOn) this.play(4); }, 1200);
},
};
Object.assign(UI, {
initMusic() {
const pop = $('#music-pop');
const open = btn => {
const r = btn.getBoundingClientRect(), w = 270;
pop.style.left = clamp(r.right - w, 8, innerWidth - w - 8) + 'px';
pop.style.top = (r.bottom + 8) + 'px';
pop.classList.add('open'); pop.setAttribute('aria-hidden', 'false');
btn.setAttribute('aria-expanded', 'true'); this._musicBtn = btn;
setTimeout(() => $('#music-toggle').focus({ preventScroll: true }), 30);
};
const close = () => { pop.classList.remove('open'); pop.setAttribute('aria-hidden', 'true'); if (this._musicBtn) this._musicBtn.setAttribute('aria-expanded', 'false'); };
this.closeMusic = close;
document.querySelectorAll('[data-music]').forEach(b => b.addEventListener('click', e => { e.stopPropagation(); pop.classList.contains('open') && this._musicBtn === b ? close() : open(b); }));
$('#music-toggle').addEventListener('click', () => Music.toggle());
$('#music-pop .mp-mix').innerHTML = mixerHTML();
$('#music-credits').addEventListener('click', () => { close(); this.openCredits(); setTimeout(() => { const t = $('#cr-musica'), sc = $('#cr-scroll'); if (t) sc.scrollTop = t.getBoundingClientRect().top - sc.getBoundingClientRect().top + sc.scrollTop - 84; }, 120); });
addEventListener('pointerdown', e => { if (pop.classList.contains('open') && !pop.contains(e.target) && !e.target.closest('[data-music]')) close(); }, true);
addEventListener('keydown', e => { if (e.key === 'Escape' && pop.classList.contains('open')) { e.stopImmediatePropagation(); close(); if (this._musicBtn) this._musicBtn.focus(); } }, true);
Music.init();
this.refreshMusic();
this.refreshAudio();
},
refreshMusic() {
if (!$('#music-pop')) return;
const on = Music.prefs.on && (Music.playing || Music.waiting);
document.querySelectorAll('[data-music]').forEach(b => { b.classList.toggle('muted', !on); b.setAttribute('aria-label', on ? 'Música: activada' : 'Música: desactivada'); b.dataset.tip = on ? 'Música: ' + MUSIC_INFO.title : 'Música desactivada'; });
const t = $('#music-toggle');
t.textContent = Music.unsupported ? 'No compatible con este navegador' : on ? (Music.playing ? 'Pausar música' : 'Comenzará al interactuar') : 'Reproducir música';
t.disabled = !!Music.unsupported; t.classList.toggle('on', on);
t.setAttribute('aria-pressed', on ? 'true' : 'false');
this.refreshAudio && this.refreshAudio();
if (Flight.on && Flight.state) Flight.refreshBar && Flight.refreshBar();
},
});
const SET_ICONS = {
general: '<circle cx="12" cy="12" r="3.2"/><path d="M12 2.8v2.6M12 18.6v2.6M2.8 12h2.6M18.6 12h2.6M5.5 5.5l1.8 1.8M16.7 16.7l1.8 1.8M5.5 18.5l1.8-1.8M16.7 7.3l1.8-1.8"/>',
audio: '<path d="M4 9.5h3.5L12 5.5v13l-4.5-4H4z"/><path d="M15.5 9a4 4 0 0 1 0 6M18 6.5a7.5 7.5 0 0 1 0 11"/>',
graphics: '<rect x="3" y="5" width="18" height="12.5" rx="2"/><path d="M8.5 21h7M12 17.5V21"/><circle cx="12" cy="11.2" r="2.6"/>',
interface: '<rect x="3" y="4" width="18" height="16" rx="2"/><path d="M3 8.5h18M8 8.5V20"/>',
controls: '<rect x="2.5" y="7" width="19" height="10.5" rx="2"/><path d="M6 10.5h1M9.5 10.5h1M13 10.5h1M16.5 10.5h1M7 14h10"/>',
simulation: '<circle cx="12" cy="12" r="3"/><ellipse cx="12" cy="12" rx="9.5" ry="4" transform="rotate(-25 12 12)"/>',
language: '<circle cx="12" cy="12" r="9"/><path d="M3 12h18M12 3c2.6 2.6 3.8 5.6 3.8 9s-1.2 6.4-3.8 9c-2.6-2.6-3.8-5.6-3.8-9S9.4 5.6 12 3z"/>',
accessibility: '<circle cx="12" cy="4.8" r="1.8"/><path d="M5 8.5l7 1.5 7-1.5M12 10v4.5l-3 6M12 14.5l3 6"/>',
};
const PCT = v => Math.round(v * 100) + ' %';
const layerGroup = (keys) => ({ get: () => keys.every(k => Settings.state.simulation.layers[k]), set: v => keys.forEach(k => { Settings.set('simulation.layers.' + k, v); UI.setLayer(k, v); }) });
const SETTINGS_SCHEMA = [
{ id: 'general', title: 'General', items: [
{ type: 'seg', path: 'general.startup', label: 'Al abrir SOLARIS', options: [['intro', 'Mostrar la pantalla de inicio'], ['explore', 'Entrar directamente al explorador']] },
{ type: 'toggle', path: 'general.confirmExit', label: 'Confirmar antes de salir de una simulación', desc: 'Pide confirmación al abandonar un vuelo en curso.' },
{ type: 'toggle', path: 'general.tips', label: 'Mostrar consejos y ayuda contextual', desc: 'Descripciones al pasar el cursor y mensajes de ayuda.' },
{ type: 'action', id: 'reset-all', label: 'Restaurar valores predeterminados', desc: 'Devuelve todos los ajustes a su configuración inicial.', danger: true },
{ type: 'note', label: 'Versión de SOLARIS: ' + SOLARIS_BUILD },
{ type: 'link', id: 'open-credits', label: 'Fuentes, créditos y licencias', desc: 'Datos de NASA, JPL y MPC; modelos 3D y música con sus atribuciones.' },
] },
{ id: 'audio', title: 'Audio', items: [
{ type: 'range', path: 'audio.master', label: 'Volumen general', desc: 'Controla todo el audio de SOLARIS.' },
{ type: 'toggle', path: 'audio.musicOn', label: 'Reproducir música' },
{ type: 'range', path: 'audio.music', label: 'Música', desc: 'Metamorphosis — Laura Platt' },
{ type: 'range', path: 'audio.ui', label: 'Efectos de interfaz', desc: 'Clics, selecciones, paneles y confirmaciones.' },
{ type: 'range', path: 'audio.nav', label: 'Navegación', desc: 'Viajes y transiciones entre objetos.' },
{ type: 'range', path: 'audio.ship', label: 'Naves y motores', desc: 'Motores, propulsores, RCS y alertas de cabina.' },
{ type: 'toggle', path: 'audio.hover', label: 'Sonidos al pasar el cursor' },
{ type: 'toggle', path: 'audio.muted', label: 'Silenciar todos los sonidos', desc: 'Tus niveles se conservan para cuando vuelvas a activar el audio.' },
] },
{ id: 'graphics', title: 'Gráficos', items: [
{ type: 'toggle', path: 'graphics.auto', label: 'Calidad gráfica automática', desc: 'Elige un perfil según tu equipo y lo reduce si la fluidez cae.' },
{ type: 'seg', path: 'graphics.preset', label: 'Perfil', preset: true, options: [['low', 'Bajo'], ['medium', 'Medio'], ['high', 'Alto'], ['ultra', 'Ultra'], ['custom', 'Personalizado']] },
{ type: 'seg', path: 'graphics.scale', label: 'Resolución de renderizado', desc: 'Más de 100 % mejora la nitidez con un coste alto de rendimiento.', options: [[0.5, '50 %'], [0.75, '75 %'], [1, '100 %'], [1.25, '125 %']] },
{ type: 'toggle', path: 'graphics.aa', label: 'Antialiasing', desc: 'Suaviza los bordes dentados.' },
{ type: 'seg', path: 'graphics.shadows', label: 'Sombras', desc: 'Eclipses de lunas y planetas y sombra de los anillos.', options: [[0, 'Desactivadas'], [1, 'Bajas'], [2, 'Medias'], [3, 'Altas']] },
{ type: 'seg', path: 'graphics.textures', label: 'Calidad de texturas', desc: 'Determina la resolución de las superficies planetarias.', options: [['low', 'Baja'], ['medium', 'Media'], ['high', 'Alta'], ['ultra', 'Ultra']] },
{ type: 'seg', path: 'graphics.distance', label: 'Distancia de renderizado', desc: 'Lunas, asteroides, cometas y naves lejanos. Los planetas siempre se dibujan.', options: [['low', 'Baja'], ['medium', 'Media'], ['high', 'Alta'], ['max', 'Máxima']] },
{ type: 'seg', path: 'graphics.lod', label: 'Nivel de detalle (LOD)', desc: 'Detalle geométrico de planetas y naves según su tamaño en pantalla.', options: [['low', 'Bajo'], ['medium', 'Medio'], ['high', 'Alto']] },
{ type: 'seg', path: 'graphics.particles', label: 'Partículas', desc: 'Estrellas, asteroides y polvo espacial.', options: [['low', 'Bajo'], ['medium', 'Medio'], ['high', 'Alto']] },
{ type: 'head', label: 'Efectos visuales' },
{ type: 'toggle', path: 'graphics.bloom', label: 'Bloom', desc: 'Resplandor alrededor de las fuentes de luz.' },
{ type: 'toggle', path: 'graphics.sunGlow', label: 'Brillo solar', desc: 'Corona y halo del Sol.' },
{ type: 'toggle', path: 'graphics.atmospheres', label: 'Atmósferas' },
{ type: 'toggle', path: 'graphics.nebulae', label: 'Vía Láctea y nebulosas' },
{ type: 'toggle', path: 'graphics.dust', label: 'Polvo espacial', desc: 'Partículas de referencia de velocidad en el modo de vuelo.' },
{ type: 'seg', path: 'graphics.sky', label: 'Calidad del fondo espacial', desc: 'Este cambio puede tardar unos segundos en aplicarse.', options: [['basic', 'Básica'], ['high', 'Alta'], ['max', 'Máxima']] },
{ type: 'seg', path: 'graphics.fps', label: 'FPS máximo', options: [[30, '30 FPS'], [60, '60 FPS'], [120, '120 FPS'], [0, 'Sin límite']] },
{ type: 'range', path: 'graphics.ambient', label: 'Luz en el lado nocturno', desc: 'Ayuda visual: en el espacio, el lado nocturno es casi negro.', min: 0, max: 0.25, step: 0.005, fmt: v => Math.round(v / 0.25 * 100) + ' %' },
] },
{ id: 'interface', title: 'Interfaz', items: [
{ type: 'seg', path: 'interface.uiScale', label: 'Escala de interfaz', options: [[0.8, '80 %'], [0.9, '90 %'], [1, '100 %'], [1.1, '110 %'], [1.25, '125 %']] },
{ type: 'range', path: 'interface.panelOpacity', label: 'Opacidad de paneles', min: 0.3, max: 0.95, step: 0.01 },
Object.assign({ type: 'toggle', id: 'labels', label: 'Etiquetas de objetos' }, layerGroup(['namesPlanets', 'namesMoons', 'namesSmall', 'namesCraft'])),
Object.assign({ type: 'toggle', id: 'orbits-ui', label: 'Mostrar órbitas', desc: 'Se recuerda para las próximas visitas.' }, layerGroup(['orbits', 'moonOrbits'])),
{ type: 'seg', path: 'interface.info', label: 'Información de las fichas', options: [['basic', 'Información básica'], ['advanced', 'Información científica avanzada']] },
{ type: 'seg', path: 'interface.hud', label: 'HUD durante el vuelo', options: [['full', 'Completo'], ['reduced', 'Reducido'], ['minimal', 'Minimalista']] },
{ type: 'seg', path: 'interface.anim', label: 'Animaciones de interfaz', desc: 'El modo reducido acorta las transiciones sin afectar ninguna función.', options: [['full', 'Activadas'], ['reduced', 'Reducidas']] },
{ type: 'toggle', path: 'interface.minimap', label: 'Minimapa' },
] },
{ id: 'controls', title: 'Controles', items: [{ type: 'keys' }] },
{ id: 'simulation', title: 'Simulación', items: [
{ type: 'seg', path: 'simulation.timeSpeed', label: 'Velocidad inicial del tiempo', options: [['paused', 'Pausado'], ['real', 'Tiempo real'], ['x10', '×10'], ['x100', '×100'], ['x1000', '×1,000']] },
{ type: 'seg', path: 'simulation.startDate', label: 'Fecha inicial', options: [['now', 'Fecha actual'], ['last', 'Última fecha utilizada']] },
Object.assign({ type: 'toggle', id: 'orbits-sim', label: 'Mostrar órbitas' }, layerGroup(['orbits', 'moonOrbits'])),
Object.assign({ type: 'toggle', id: 'moons', label: 'Mostrar lunas' }, layerGroup(['moons'])),
Object.assign({ type: 'toggle', id: 'asteroids', label: 'Mostrar asteroides', desc: 'Cinturón principal y cinturón de Kuiper.' }, layerGroup(['asteroids', 'kuiper'])),
Object.assign({ type: 'toggle', id: 'craft', label: 'Mostrar satélites artificiales' }, layerGroup(['craft'])),
Object.assign({ type: 'toggle', id: 'const', label: 'Constelaciones' }, layerGroup(['constLines', 'constNames'])),
{ type: 'seg', path: 'simulation.scale', label: 'Escala', options: [['visual', 'Escala visual'], ['real', 'Distancias más realistas'], ['edu', 'Educativa']] },
{ type: 'toggle', path: 'simulation.rotCap', label: 'Ralentizar la rotación a velocidades altas', desc: 'Evita el parpadeo cuando un día dura menos de un fotograma.' },
] },
{ id: 'accessibility', title: 'Accesibilidad', items: [
{ type: 'toggle', path: 'accessibility.reducedMotion', label: 'Reducir movimiento', desc: 'Viajes de cámara breves, sin vibración ni efectos de aceleración.' },
{ type: 'toggle', path: 'accessibility.contrast', label: 'Alto contraste' },
{ type: 'seg', path: 'accessibility.text', label: 'Tamaño de texto', options: [['small', 'Pequeño'], ['normal', 'Normal'], ['large', 'Grande']] },
{ type: 'toggle', path: 'accessibility.reduceFlashes', label: 'Reducir efectos luminosos', desc: 'Atenúa el bloom, el brillo solar y el resplandor de los motores.' },
{ type: 'toggle', path: 'audio.hover', label: 'Sonidos al pasar el cursor' },
] },
];
Object.assign(UI, {
setCat: 'general',
openSettings(cat) {
if (cat) this.setCat = cat;
const m = $('#settings'); this._setFrom = document.activeElement;
this.renderSettings();
m.classList.add('open'); m.setAttribute('aria-hidden', 'false');
setTimeout(() => { const b = m.querySelector('.set-nav .on'); if (b) b.focus({ preventScroll: true }); }, 40);
},
closeSettings() {
const m = $('#settings'); if (!m.classList.contains('open')) return;
this.cancelRebind();
m.classList.remove('open'); m.setAttribute('aria-hidden', 'true');
if (this._setFrom && this._setFrom.focus) this._setFrom.focus({ preventScroll: true });
},
itemValue(it) { return it.get ? it.get() : it.path ? Settings.get(it.path) : undefined; },
renderSettings() {
const cat = SETTINGS_SCHEMA.find(c => c.id === this.setCat) || SETTINGS_SCHEMA[0];
$('#set-nav').innerHTML = SETTINGS_SCHEMA.map(c => `<button role="tab" aria-selected="${c.id === cat.id}" class="${c.id === cat.id ? 'on' : ''}" data-cat="${c.id}"><svg viewBox="0 0 24 24" aria-hidden="true">${SET_ICONS[c.id]}</svg><span>${t(c.title)}</span></button>`).join('');
const g = Settings.state.graphics;
$('#set-body').innerHTML = `<h2 id="set-cat-title">${t(cat.title)}</h2>` + cat.items.map((it, i) => {
const id = 'set-' + cat.id + '-' + i, v = this.itemValue(it), desc = it.desc ? `<p class="set-desc">${t(it.desc)}</p>` : '';
if (it.type === 'head') return `<h3 class="set-head">${t(it.label)}</h3>`;
if (it.type === 'note') return `<p class="set-note">${t(it.label)}</p>`;
if (it.type === 'keys') return this.keysHTML();
if (it.type === 'action' || it.type === 'link') return `<div class="set-row"><div class="set-lbl"><b>${t(it.label)}</b>${desc}</div><button class="${it.danger ? 'set-btn danger' : 'set-btn'}" data-act="${it.id}">${t(it.type === 'link' ? 'Abrir' : 'Restaurar')}</button></div>`;
if (it.type === 'toggle') return `<label class="set-row set-tog" for="${id}"><div class="set-lbl"><b>${t(it.label)}</b>${desc}</div><input type="checkbox" id="${id}" data-i="${i}" ${v ? 'checked' : ''}><span class="sw" aria-hidden="true"></span></label>`;
if (it.type === 'range') { const mn = it.min == null ? 0 : it.min, mx = it.max == null ? 1 : it.max, st = it.step || 0.01, f = it.fmt || PCT;
return `<div class="set-row set-range"><div class="set-lbl"><label for="${id}"><b>${t(it.label)}</b></label>${desc}</div><div class="set-ctl"><input type="range" id="${id}" data-i="${i}" min="${mn}" max="${mx}" step="${st}" value="${v}"><output>${f(v)}</output></div></div>`; }
if (it.type === 'seg') {
const warn = it.preset && g.preset === 'ultra' ? `<p class="set-warn">${t('Puede afectar el rendimiento en equipos menos potentes.')}</p>` : '';
return `<div class="set-row set-seg"><div class="set-lbl"><b id="${id}-l">${t(it.label)}</b>${desc}${warn}</div><div class="seg" role="radiogroup" aria-labelledby="${id}-l">${it.options.map(([ov, ol]) =>
`<button role="radio" aria-checked="${ov === v}" class="${ov === v ? 'on' : ''}" data-i="${i}" data-v='${JSON.stringify(ov)}' ${it.preset && ov === 'custom' ? 'disabled' : ''}>${t(ol)}</button>`).join('')}</div></div>`;
}
return '';
}).join('');
this.renderKeyHelp();
},
keysHTML() {
const row = (ctx, [a]) => { const k = Keys.code(ctx, a); return `<div class="key-row"><span>${t(Keys.describe(ctx, a))}</span><button class="kbd-btn" data-ctx="${ctx}" data-act-key="${a}" aria-label="${t('Cambiar la tecla de «{0}»', t(Keys.describe(ctx, a)))}">${Keys.label(k)}</button></div>`; };
const fixed = (rows) => rows.map(([l, k]) => `<div class="key-row fixed"><span>${t(l)}</span><kbd>${t(k)}</kbd></div>`).join('');
return `<p class="set-desc">${t('Haz clic en una tecla y pulsa la nueva. Si ya estaba en uso, se intercambia para que ninguna acción quede sin control.')}</p>
<p id="key-msg" class="set-note" role="status" aria-live="polite"></p>
<h3 class="set-head">${t('Exploración')}</h3><div class="key-grid">${fixed([['Rotar cámara', 'Arrastrar'], ['Zoom', 'Rueda'], ['Seleccionar objeto', 'Clic'], ['Sol, planetas y Plutón', '0–9']])}${KEY_ACTIONS.explore.map(r => row('explore', r)).join('')}</div>
<h3 class="set-head">${t('Nave espacial')}</h3><div class="key-grid">${fixed([['Empuje 25, 50, 75 y 100 %', '1–4'], ['Cámaras', 'Mayús + 1–6'], ['Pausa y ayuda', 'Esc']])}${KEY_ACTIONS.flight.map(r => row('flight', r)).join('')}</div>
<div class="set-row"><div class="set-lbl"><b>${t('Restaurar controles predeterminados')}</b></div><button class="set-btn" data-act="reset-keys">${t('Restaurar')}</button></div>`;
},
initSettings() {
const m = $('#settings');
m.addEventListener('click', e => {
if (e.target === m) return this.closeSettings();
const nav = e.target.closest('[data-cat]'); if (nav) { this.setCat = nav.dataset.cat; this.renderSettings(); $('#set-body').scrollTop = 0; return; }
const seg = e.target.closest('.seg [data-i]');
if (seg && !seg.disabled) {
const it = this.curItems()[+seg.dataset.i], v = JSON.parse(seg.dataset.v);
if (it.preset) { Settings.applyPreset(v); Settings.state.graphics.auto = false; Settings.save(); }
else if (it.set) it.set(v); else Settings.set(it.path, v);
this.renderSettings(); return;
}
const kb = e.target.closest('[data-act-key]'); if (kb) return this.startRebind(kb);
const act = e.target.closest('[data-act]'); if (!act) return;
if (act.dataset.act === 'reset-all') this.confirmDialog(t('¿Restaurar todos los ajustes?'), t('Se borrarán tus preferencias de audio, gráficos, interfaz, controles e idioma.'), t('Restaurar'), t('Cancelar')).then(ok => { if (ok) { Settings.reset(); this.renderSettings(); this.toast(t('Ajustes restaurados')); } });
if (act.dataset.act === 'reset-keys') { Settings.resetControls(); this.renderSettings(); this.toast(t('Controles restaurados')); }
if (act.dataset.act === 'open-credits') { this.closeSettings(); this.openCredits(); }
});
m.addEventListener('change', e => { const el = e.target; if (!el.matches('input[type="checkbox"][data-i]')) return; const it = this.curItems()[+el.dataset.i]; if (it.set) it.set(el.checked); else Settings.set(it.path, el.checked); if (it.path === 'graphics.auto') this.renderSettings(); });
m.addEventListener('input', e => { const el = e.target; if (!el.matches('input[type="range"][data-i]')) return; const it = this.curItems()[+el.dataset.i], v = +el.value; Settings.set(it.path, v); el.nextElementSibling.textContent = (it.fmt || PCT)(v); });
$('#set-close').addEventListener('click', () => this.closeSettings());
$('#set-done').addEventListener('click', () => this.closeSettings());
$('#set-reset').addEventListener('click', () => this.confirmDialog(t('¿Restaurar todos los ajustes?'), t('Se borrarán tus preferencias de audio, gráficos, interfaz, controles e idioma.'), t('Restaurar'), t('Cancelar')).then(ok => { if (ok) { Settings.reset(); this.renderSettings(); this.toast(t('Ajustes restaurados')); } }));
document.querySelectorAll('[data-open-settings]').forEach(b => b.addEventListener('click', () => this.openSettings()));
addEventListener('keydown', e => { if (e.key === 'Escape' && m.classList.contains('open') && !this.rebind) { e.stopImmediatePropagation(); this.closeSettings(); } }, true);
Settings.on(() => { if (m.classList.contains('open') && !this._rendering) { this._rendering = true; requestAnimationFrame(() => { this._rendering = false; if (!this.rebind && !document.activeElement.matches('input[type="range"]')) this.renderSettings(); }); } });
document.addEventListener('solaris:lang', () => this.onLanguage());
this.renderKeyHelp();
},
curItems() { return (SETTINGS_SCHEMA.find(c => c.id === this.setCat) || SETTINGS_SCHEMA[0]).items; },
startRebind(btn) {
this.cancelRebind();
this.rebind = { btn, ctx: btn.dataset.ctx, act: btn.dataset.actKey, old: btn.textContent };
btn.textContent = t('Pulsa una tecla…'); btn.classList.add('listening');
this.rebindH = e => {
e.preventDefault(); e.stopImmediatePropagation();
if (e.key === 'Escape') return this.cancelRebind();
const key = Keys.norm(e), r = Keys.assign(this.rebind.ctx, this.rebind.act, key);
const msg = $('#key-msg');
if (!r.ok) { if (msg) msg.textContent = t('La tecla «{0}» está reservada para la navegación básica.', Keys.label(key)); return; }
const swapped = r.swapped ? t('Se intercambió con «{0}».', t(Keys.describe(this.rebind.ctx, r.swapped))) : '';
this.cancelRebind(true); this.renderSettings();
const m2 = $('#key-msg'); if (m2) m2.textContent = t('Asignada: {0}.', Keys.label(key)) + (swapped ? ' ' + swapped : '');
};
addEventListener('keydown', this.rebindH, true);
},
cancelRebind(done) {
if (!this.rebind) return;
removeEventListener('keydown', this.rebindH, true);
if (!done) { this.rebind.btn.textContent = this.rebind.old; this.rebind.btn.classList.remove('listening'); }
this.rebind = null;
},
renderKeyHelp() {
const K = (c, a) => `<kbd>${Keys.label(Keys.code(c, a))}</kbd>`;
const dl = rows => '<dl>' + rows.map(([k, d]) => `<dt>${k}</dt><dd>${t(d)}</dd>`).join('') + '</dl>';
const ex = a => [K('explore', a), Keys.describe('explore', a)];
const fl = a => [K('flight', a), Keys.describe('flight', a)];
const help = $('#help-keys');
if (help) help.innerHTML = dl([['<kbd>' + t('Arrastrar') + '</kbd>', 'Girar la cámara'], ['<kbd>' + t('Rueda') + '</kbd>', 'Acercar y alejar'], ['Doble clic', 'Viajar muy cerca'], ex('rotL'), ex('zoomIn'), ex('zoomOut')])
+ dl([['<kbd>0</kbd>–<kbd>9</kbd>', 'Sol, planetas y Plutón'], ex('search'), ex('explore'), ex('center'), ex('follow'), ex('overview'), ex('sun')])
+ dl([ex('pause'), ex('faster'), ex('slower'), ex('now'), ex('names'), ex('orbits'), ex('constel'), ex('minimap'), ex('cine'), ex('tours'), ['<kbd>Esc</kbd>', 'Cerrar o salir'], ['<kbd>' + t('Inicio') + '</kbd>', 'Pantalla inicial'], ex('pilot')]);
const ge = $('#guide-keys-explore'), gf = $('#guide-keys-flight');
if (ge) ge.innerHTML = dl([['<kbd>0</kbd>–<kbd>9</kbd>', 'Sol, planetas y Plutón'], ex('search'), ex('explore'), ex('center'), ex('follow'), ex('overview'), ex('sun'), ex('zoomIn'), ex('zoomOut')])
+ dl([ex('pause'), ex('faster'), ex('slower'), ex('now'), ex('names'), ex('orbits'), ex('constel'), ex('minimap'), ex('cine'), ex('tours'), ex('pilot'), ['<kbd>?</kbd>', 'Atajos rápidos'], ['<kbd>Esc</kbd>', 'Cerrar o salir']]);
if (gf) gf.innerHTML = dl([fl('thrUp'), fl('thrDown'), ['<kbd>1</kbd>–<kbd>4</kbd>', 'Empuje 25, 50, 75 y 100 %'], fl('thr0'), fl('brake'), fl('pitchUp'), fl('yawL'), fl('rollL'), fl('left'), fl('up')])
+ dl([fl('target'), fl('intercept'), fl('aim'), fl('match'), fl('hold'), fl('orbit'), fl('autobrake'), fl('stab'), fl('mode'), fl('cam'), fl('warpUp'), fl('sound'), ['<kbd>Esc</kbd>', 'Pausa y ayuda']]);
},
confirmDialog(title, text, ok, cancel) {
return new Promise(res => {
const d = $('#confirm'); $('#cf-title').textContent = title; $('#cf-text').textContent = text; $('#cf-ok').textContent = ok; $('#cf-cancel').textContent = cancel;
d.classList.add('open'); d.setAttribute('aria-hidden', 'false'); const prev = document.activeElement;
const done = v => { d.classList.remove('open'); d.setAttribute('aria-hidden', 'true'); $('#cf-ok').onclick = $('#cf-cancel').onclick = null; removeEventListener('keydown', kh, true); if (prev && prev.focus) prev.focus({ preventScroll: true }); res(v); };
const kh = e => { if (e.key === 'Escape') { e.preventDefault(); e.stopImmediatePropagation(); done(false); } };
$('#cf-ok').onclick = () => done(true); $('#cf-cancel').onclick = () => done(false);
addEventListener('keydown', kh, true);
setTimeout(() => $('#cf-cancel').focus({ preventScroll: true }), 30);
});
},
onLanguage() {
this.relabel();
if ($('#settings').classList.contains('open')) this.renderSettings();
this.renderKeyHelp();
if (this.sel) this.renderInfo(this.sel);
if (this.section) this.openSection(this.section);
if (this.tour && !this.tour.done) { const st = this.tour.def.stops[this.tour.i]; this.renderTourCard(st, World.byId[st.id]); this.refreshTourBar(); }
if ($('#tours-list')) this.renderTourCards && this.renderTourCards();
this.refreshTime(); this.tickTime && this.tickTime(); this.refreshMusic(); this.updateScaleBadge && this.updateScaleBadge();
if (Flight.on) { if (Flight.state === 'hangar') Flight.renderHangar(); else Flight.refreshBar(); }
},
relabel() {
for (const rb of World.rb) { const L = this.labels[rb.id]; if (L && L.n) { L.n.textContent = rb.def.short || rb.def.name; L.w = 0; } }
for (const c of this.constLabels) { const src = CONSTELLATIONS.find(x => x.id === c.c.id); if (src) { c.c.name = src.name; c.el.textContent = src.name; } }
},
});
const SettingsApply = {
timers: {},
debounce(k, ms, f) { clearTimeout(this.timers[k]); this.timers[k] = setTimeout(f, ms); },
boot() {
const st = Settings.state || Settings.load();
if (st.graphics.auto) { const p = Settings.detectQuality(); Object.assign(st.graphics, GFX_PRESETS[p]); st.graphics.preset = p; }
this.syncRuntime();
this.applyCSS();
GLX.setTexQuality(st.graphics.textures);
Settings.on((path, v) => this.onChange(path, v));
addEventListener('pagehide', () => { this.saveDate(); Settings.flush(); });
document.addEventListener('visibilitychange', () => { if (document.hidden) Settings.flush(); });
},
syncRuntime() {
const st = Settings.state, sim = st.simulation;
S.scale = sim.scale; ScaleState.set(sim.scale, false);
S.ambient = st.graphics.ambient; S.rotCap = sim.rotCap; S.hzModel = HZ_MODELS[sim.hzModel] ? sim.hzModel : 'conservador';
Object.assign(S.layers, sim.layers);
UI.mmOn = st.interface.minimap && innerWidth >= 760;
},
applyCSS() {
const st = Settings.state, h = document.documentElement, I = st.interface, A = st.accessibility;
h.classList.toggle('anim-reduced', I.anim === 'reduced' || A.reducedMotion);
h.classList.toggle('hc', A.contrast);
h.classList.remove('ts-small', 'ts-large'); if (A.text !== 'normal') h.classList.add('ts-' + A.text);
h.style.setProperty('--uiz', I.uiScale);
const a = I.panelOpacity;
h.style.setProperty('--glass', `rgba(10, 15, 27, ${a})`);
h.style.setProperty('--glass-2', `rgba(12, 18, 32, ${Math.min(0.97, a + 0.24)})`);
const app = document.getElementById('app');
if (app) { app.classList.toggle('hud-reduced', I.hud === 'reduced'); app.classList.toggle('hud-min', I.hud === 'minimal'); }
},
onChange(path, v) {
const st = Settings.state;
if (path === '*') { this.syncRuntime(); this.applyCSS(); App.applyQuality(true); GLX.setTexQuality(st.graphics.textures); this.applyAtmos(); this.rebakeSky(); SFX.apply(); Music.applyLevel(); UI.toggleMinimap(UI.mmOn); I18N.set(st.language); UI.refreshAudio(); this.refreshViews(); return; }
const [sec, key, sub] = path.split('.');
if (sec === 'audio') {
SFX.apply();
if (key === 'music') Music.applyLevel();
if (key === 'musicOn') { if (v) Music.play(2.2); else Music.pause(); }
UI.refreshAudio(); UI.refreshMusic();
} else if (sec === 'graphics') {
if (['scale', 'aa', 'preset'].includes(key)) this.debounce('q', 250, () => App.applyQuality(true));
if (key === 'textures' || key === 'preset') GLX.setTexQuality(st.graphics.textures);
if (key === 'sky' || key === 'preset') this.debounce('sky', 400, () => this.rebakeSky());
if (key === 'atmospheres' || key === 'preset') this.applyAtmos();
if (key === 'ambient') S.ambient = v;
if (key === 'auto' && v) { const p = Settings.detectQuality(); Settings.applyPreset(p); st.graphics.auto = true; Settings.save(); }
} else if (sec === 'interface') {
this.applyCSS();
if (key === 'minimap') UI.toggleMinimap(v);
if (key === 'info' && UI.sel) UI.renderInfo(UI.sel);
} else if (sec === 'accessibility') {
this.applyCSS();
} else if (sec === 'simulation') {
if (key === 'layers') { S.layers[sub] = v; document.querySelectorAll(`[data-layer="${sub}"]`).forEach(el => { el.checked = v; }); }
if (key === 'scale') UI.setScale(v);
if (key === 'rotCap') S.rotCap = v;
} else if (sec === 'language') {
I18N.set(v);
}
if (sec === 'controls' || sec === 'language') UI.renderKeyHelp && UI.renderKeyHelp();
},
applyAtmos() {
const on = Settings.state.graphics.atmospheres;
(World.rb || []).forEach(rb => {
const A = rb.def.vis && rb.def.vis.atm; if (!A || !rb.u) return;
rb.u.u_atm = on ? [...A.color.map(x => x * 0.9), A.k] : [0, 0, 0, 0];
});
},
rebakeSky() {
if (!World.sky) return;
const size = { basic: 512, high: 1024, max: 2048 }[Settings.state.graphics.sky];
if (size === World.skySize) return;
UI.toast(I18N.t('Este cambio puede tardar unos segundos en aplicarse.'), 2200);
setTimeout(() => { World.bakeSky(); GLX.resize(innerWidth, innerHeight, App.dpr, true); }, 60);
},
refreshViews() {
if (UI.sel) UI.renderInfo(UI.sel);
if (UI.section) UI.openSection(UI.section);
UI.refreshTime && UI.refreshTime();
},
autoCheck(dt) {
const g = Settings.state.graphics; if (!g.auto || document.hidden || !UI.started) return;
this.acc = this.acc || { t: 0, n: 0, s: 0, cool: 6 };
const a = this.acc; a.cool -= dt; if (a.cool > 0) return;
a.t += dt; a.n++; if (a.t < 5) return;
const fps = a.n / a.t; a.t = 0; a.n = 0;
const order = ['low', 'medium', 'high'], i = order.indexOf(g.preset);
if (fps < 26 && i > 0) { Settings.applyPreset(order[i - 1]); g.auto = true; Settings.save(); a.cool = 12; UI.hint(I18N.t('Calidad gráfica ajustada automáticamente a «{0}» para mantener la fluidez.', I18N.t(GFX_NAMES[order[i - 1]]))); }
},
saveDate() { if (!Flight.on && window.Time) { Settings.state.simulation.lastJD = Time.jd; Settings.save(); } },
};
const GFX_NAMES = { low: 'Bajo', medium: 'Medio', high: 'Alto', ultra: 'Ultra', custom: 'Personalizado' };
const CINE_ICON = {
prev: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M15 6l-6 6 6 6"/></svg>',
next: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M9 6l6 6-6 6"/></svg>',
pause: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M8.5 6v12M15.5 6v12"/></svg>',
play: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M8 5.5v13l10.5-6.5z"/></svg>',
ui: '<svg viewBox="0 0 24 24" aria-hidden="true"><rect x="3" y="5" width="18" height="14" rx="2"/><path d="M3 9.5h18M8 9.5V19"/></svg>',
exit: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M6 6l12 12M18 6L6 18"/></svg>',
film: '<svg viewBox="0 0 24 24" aria-hidden="true"><rect x="3" y="5" width="18" height="14" rx="2"/><path d="M7 5v14M17 5v14M3 9.5h4M3 14.5h4M17 9.5h4M17 14.5h4"/></svg>',
};
Object.assign(UI, {
tourCine: false,
initTourCine() {
const c = $('#cine-ctl');
c.innerHTML = `
<button data-c="prev" aria-label="Parada anterior" data-tip="Anterior (←)">${CINE_ICON.prev}</button>
<button data-c="play" class="cc-play" aria-label="Pausar recorrido" data-tip="Pausar (espacio)">${CINE_ICON.pause}</button>
<button data-c="next" aria-label="Parada siguiente" data-tip="Siguiente (→)">${CINE_ICON.next}</button>
<span class="cc-sep"></span>
<span class="cc-step" id="cc-step"></span>
<span class="cc-track"><i id="cc-prog"></i></span>
<span class="cc-sep"></span>
<button data-c="ui" aria-label="Salir del Modo Cine y seguir el recorrido" data-tip="Salir del Modo Cine">${CINE_ICON.ui}</button>
<button data-c="exit" aria-label="Salir del recorrido" data-tip="Salir del recorrido (Esc)">${CINE_ICON.exit}</button>`;
c.addEventListener('click', e => {
const b = e.target.closest('[data-c]'); if (!b || !this.tour) return;
const a = b.dataset.c;
if (a === 'prev') this.stepTour(-1); else if (a === 'next') this.stepTour(1);
else if (a === 'play') this.pauseTour(!this.tour.paused, true);
else if (a === 'ui') this.exitTourCine();
else if (a === 'exit') this.stopTour();
this.cineWake(); $('#gl').focus({ preventScroll: true });
});
['pointermove', 'pointerdown', 'keydown', 'touchstart', 'wheel'].forEach(t => addEventListener(t, () => { if (this.tourCine) this.cineWake(); }, { passive: true, capture: true }));
$('#tb-cine').addEventListener('click', () => { this.enterTourCine(); $('#gl').focus({ preventScroll: true }); });
},
enterTourCine() {
if (this.tourCine || !this.tour) return;
this.tourCine = true;
this.openSection(null); this.toggleHelp(false);
if (this.sel) this.select(null, { keepCam: true, tour: true });
$('#app').classList.add('tour-cine');
SFX.apply();
this.cineWake();
this.refreshCine();
const T = this.tour;
if (T && T.phase !== 'intro' && !T.done && this.tourRB) setTimeout(() => this.cineCaption(T.def.stops[T.i], this.tourRB), 900);
},
exitTourCine() {
if (!this.tourCine) return;
this.tourCine = false;
const app = $('#app'); app.classList.add('cine-out'); app.classList.remove('tour-cine', 'cine-active');
clearTimeout(this.outT); this.outT = setTimeout(() => app.classList.remove('cine-out'), 1300);
$('#cine-cap').classList.remove('on'); clearTimeout(this.capT); clearTimeout(this.wakeT);
SFX.apply();
if (this.tour && this.tour.phase === 'hold') Cam.auto = 0.05;
},
cineWake() {
$('#app').classList.add('cine-active');
clearTimeout(this.wakeT);
this.wakeT = setTimeout(() => { if (this.tourCine && this.tour && !this.tour.paused) $('#app').classList.remove('cine-active'); }, 3200);
},
refreshCine() {
const T = this.tour; if (!T || !this.tourCine) return;
const n = T.def.stops.length;
$('#cc-step').textContent = (T.i + 1) + ' / ' + n;
const p = $('#cine-ctl .cc-play');
p.innerHTML = T.paused ? CINE_ICON.play : CINE_ICON.pause;
p.setAttribute('aria-label', T.paused ? 'Reanudar recorrido' : 'Pausar recorrido');
if (T.paused) this.cineWake();
},
cineCaption(stop, rb) {
if (!this.tourCine || !rb) return;
const d = rb.def, cap = $('#cine-cap');
let title = (stop.kind === 'system' && stop.group) ? stop.group : d.name, sub = '';
const F = rb.isCraft ? null : Info.fields(rb), all = F ? [...F.phys, ...F.orb] : [];
const get = re => { const f = all.find(x => re.test(x.label)); return f ? f.val : null; };
if (rb.isCraft) sub = d.sub + (d.craft && d.craft.agency ? ' • ' + d.craft.agency : '');
else if (stop.kind === 'system') { const m = get(/^Número de lunas/); const n = m ? m.replace(/\s*\(.*\)/, '') : null; sub = n ? n + (n === '1' ? ' luna conocida' : ' lunas conocidas') : TYPE_LABEL[d.type]; }
else if (rb.isSun) sub = 'Estrella • el centro del Sistema Solar';
else if (d.type === 'moon' && rb.parent) { const dm = get(/^Distancia media a/); sub = 'Luna de ' + rb.parent.def.name + (dm ? ' • ' + dm + ' de ' + rb.parent.def.name : ''); }
else {
const r = Info.live(rb).find(x => x.label === 'Distancia actual al Sol');
sub = TYPE_LABEL[d.type] + (r && r.km ? ' • ' + fmtKm(r.km) + ' del Sol' : '');
}
cap.innerHTML = `<b>${esc(title.toUpperCase())}</b><span>${esc(sub)}</span>`;
cap.classList.add('on');
clearTimeout(this.capT);
this.capT = setTimeout(() => cap.classList.remove('on'), 6500);
},
hideCineCaption() { clearTimeout(this.capT); $('#cine-cap').classList.remove('on'); },
});
const JD_UNIX = 2440587.5;
const jdToDate = jd => new Date((jd - JD_UNIX) * 86400000);
const dateToJD = d => d.getTime() / 86400000 + JD_UNIX;
ICON.herramientas = '<svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M4 20h16M7 16V9M12 16V5M17 16v-4"/></svg>';
UI.SECTIONS.splice(UI.SECTIONS.findIndex(s => s[0] === 'capas'), 0, ['herramientas', 'Herramientas']);
Object.assign(UI, {
setSimDate(jd, opts) {
opts = opts || {};
Time.jd = jd; Time.live = false;
if (opts.pause) Time.paused = true;
this.refreshTime();
if (opts.toast !== false) {
const y = jdToDate(jd).getUTCFullYear();
this.toast(y < 1800 || y > 2050 ? 'Fuera de 1800–2050 la precisión de las posiciones disminuye.' : 'Posiciones aproximadas para el ' + this.dateStr(jd));
}
TL.refresh();
},
isSimulatedDate() { return !Time.live || Math.abs(Time.jd - dateToJD(new Date())) > 0.5; },
});
const _sectionHTML = UI.sectionHTML;
UI.sectionHTML = function (id) {
if (id !== 'herramientas') return _sectionHTML.call(this, id);
const row = (t, ico, title, sub) => `<button class="cta-row" data-tool="${t}">${ico}<span><b>${title}</b><small>${sub}</small></span></button>`;
return `<p class="lead">Herramientas para comparar, medir y viajar en el tiempo. Todas usan la misma fecha y los mismos datos de la simulación.</p>
<div class="list">${row('compare', ICON.herramientas, 'Comparar cuerpos celestes', 'Hasta 4 objetos: datos, tamaño real y gravedad')}
${row('timeline', ICON.calendar, 'Línea del tiempo', 'Recorre fechas y eventos de la exploración espacial')}
${(UI.TOOL_ROWS || []).map(r => row(...r)).join('')}</div>`;
};
document.addEventListener('click', e => {
const t = e.target.closest('[data-tool]'); if (!t) return;
({ compare: () => Compare.open(), timeline: () => TL.toggle(true) })[t.dataset.tool]?.() ?? (UI.TOOL_ACTIONS || {})[t.dataset.tool]?.();
});
const Compare = {
list: [], view: 'tabla', MAX: 4,
can(rb) { return rb && !rb.isCraft && ['star', 'planet', 'dwarf', 'tno', 'moon', 'asteroid', 'comet'].includes(rb.def.type); },
has(id) { return this.list.includes(id); },
add(id) {
const rb = World.byId[id]; if (!this.can(rb) || this.has(id)) return;
if (this.list.length >= this.MAX) { UI.toast('El comparador admite hasta ' + this.MAX + ' objetos.'); return; }
this.list.push(id); this.changed();
if (this.list.length === 1) UI.toast(rb.def.name + ' añadido. Añade otro objeto desde su ficha o desde el comparador.');
},
remove(id) { this.list = this.list.filter(x => x !== id); this.changed(); },
clear() { this.list = []; this.changed(); },
changed() {
this.paintTray();
if (UI.sel) { const b = $('#info-cmp'); if (b) b.outerHTML = this.pillHTML(UI.sel); }
if (this.isOpen()) this.render();
},
pillHTML(rb) {
if (!this.can(rb)) return '';
const on = this.has(rb.id);
return `<button class="pill${on ? ' on' : ''}" id="info-cmp" data-cmp="${rb.id}" aria-pressed="${on}">${on ? '✓ En comparación' : '+ Comparar'}</button>`;
},
paintTray() {
let t = $('#cmp-tray');
if (!t) { t = document.createElement('div'); t.id = 'cmp-tray'; t.className = 'ui'; $('#app').appendChild(t); }
const n = this.list.length; t.classList.toggle('on', n > 0 && !this.isOpen());
t.innerHTML = n ? `<span class="ct-dots">${this.list.map(id => `<i class="dot" style="--c:${BODY[id].color}"></i>`).join('')}</span><span>${this.list.map(id => esc(BODY[id].short || BODY[id].name)).join(' · ')}</span>
<button class="cta" data-cmp-open>${n >= 2 ? 'Comparar (' + n + ')' : 'Añade otro objeto'}</button><button class="icon-btn" data-cmp-clear aria-label="Vaciar comparación">&times;</button>` : '';
},
isOpen() { const p = $('#cmp'); return !!p && p.classList.contains('open'); },
open() {
let p = $('#cmp');
if (!p) {
p = document.createElement('section'); p.id = 'cmp'; p.className = 'tool-panel ui'; p.setAttribute('role', 'dialog'); p.setAttribute('aria-labelledby', 'cmp-title');
$('#app').appendChild(p);
p.addEventListener('click', e => {
const a = e.target.closest('[data-c]'); if (!a) return; const v = a.dataset.c, id = a.dataset.id;
if (v === 'close') this.close(); else if (v === 'remove') this.remove(id);
else if (v === 'locate') { this.close(); UI.select(World.byId[id], { fly: true }); }
else if (v === 'view') { this.view = a.dataset.v; this.render(); }
else if (v === 'add') { this.add(id); $('#cmp-q').value = ''; }
});
p.addEventListener('input', e => { if (e.target.id === 'cmp-q') this.renderHits(); });
addEventListener('keydown', e => { if (e.key === 'Escape' && this.isOpen()) { e.stopImmediatePropagation(); this.close(); } }, true);
}
UI.openSection(null);
p.classList.add('open'); this.render(); this.paintTray();
setTimeout(() => { const q = $('#cmp-q'); if (q && this.list.length < 2) q.focus(); }, 60);
},
close() { const p = $('#cmp'); if (p) p.classList.remove('open'); this.paintTray(); },
ROWS: [
['Tipo', null], ['Diámetro', /^Diámetro medio/, 'diam'], ['Radio', /^Radio medio/], ['Masa', /^Masa$/, 'mass'],
['Gravedad superficial', /^Gravedad superficial/, 'grav'], ['Temperatura', /^Temperatura aproximada/], ['Duración del día', /^Duración de un día/],
['Periodo orbital', /^(Duración de un año|Periodo orbital)/, 'per'], ['Distancia media al Sol', /^Distancia media al Sol/, 'au'],
['Velocidad orbital', /^Velocidad orbital media/, 'v'], ['Número de lunas', /^Número de lunas/], ['Densidad', /^Densidad media/, 'dens'],
['Composición', /^Composición/], ['Atmósfera', /^Tipo de atmósfera/],
],
num(rb, k) {
const d = rb.def, i = d.info;
if (k === 'diam') return d.R ? 2 * d.R : null;
if (k === 'mass') return d.mass || null;
if (k === 'grav') return this.grav(rb);
if (k === 'dens') return d.mass && d.R ? d.mass * 1000 / (4 / 3 * Math.PI * Math.pow(d.R * 1e5, 3)) : null;
if (k === 'per') { const P = Info.periodDays(rb); return P ? P.v : null; }
if (k === 'v') return typeof i.v === 'number' ? i.v : null;
if (k === 'au') { const f = this.field(rb, /^Distancia media al Sol/); const m = f && String(f.val).match(/([\d.,]+)\s*UA/); return m ? parseFloat(m[1].replace(/,/g, '')) : null; }
return null;
},
grav(rb) { const d = rb.def; return d.info.grav != null ? d.info.grav : d.mass && d.R ? G_CONST * d.mass / Math.pow(d.R * 1000, 2) : null; },
field(rb, re) { const F = rb._cmpF || (rb._cmpF = Info.fields(rb)); return [...F.phys, ...F.orb, ...F.comp].find(f => re.test(f.label) && f.val != null); },
bars(vals) {
const ok = vals.filter(v => v != null && v > 0); if (ok.length < 2) return null;
const mx = Math.max(...ok), mn = Math.min(...ok), log = mx / mn > 40;
return { log, k: v => v == null || v <= 0 ? 0 : log ? Math.max(0.03, (Math.log10(v) - Math.log10(mn) + 0.3) / (Math.log10(mx) - Math.log10(mn) + 0.3)) : Math.max(0.03, v / mx) };
},
render() {
const p = $('#cmp'); if (!p) return;
const rbs = this.list.map(id => World.byId[id]).filter(Boolean); rbs.forEach(rb => { rb._cmpF = null; });
const tabs = [['tabla', 'Datos'], ['tamano', 'Comparar tamaño'], ['gravedad', 'Gravedad']];
p.innerHTML = `
<header class="tp-head"><div><span class="tp-kicker">Herramienta</span><h2 id="cmp-title">Comparar cuerpos celestes</h2></div><button class="icon-btn" data-c="close" aria-label="Cerrar">&times;</button></header>
<div class="cmp-cards">${rbs.map(rb => `<div class="cmp-card"><i class="dot" style="--c:${rb.def.color}"></i><b>${esc(rb.def.name)}</b><small>${esc(TYPE_LABEL[rb.def.type])}${rb.def.type === 'moon' && rb.parent ? ' de ' + esc(rb.parent.def.name) : ''}</small>
<span class="cmp-card-act"><button class="txt-btn" data-c="locate" data-id="${rb.id}">Localizar</button><button class="txt-btn" data-c="remove" data-id="${rb.id}" aria-label="Quitar ${esc(rb.def.name)}">Quitar</button></span></div>`).join('')}
${rbs.length < this.MAX ? `<div class="cmp-card cmp-add"><input id="cmp-q" type="search" placeholder="Añadir un objeto…" autocomplete="off" aria-label="Buscar un objeto para comparar"><ul id="cmp-hits" role="listbox"></ul></div>` : ''}</div>
${rbs.length < 2 ? `<p class="note cmp-empty">Elige al menos dos objetos para compararlos. Puedes añadirlos aquí o con «+ Comparar» en la ficha de cada objeto.</p>` : `
<div class="seg tp-tabs" role="tablist">${tabs.map(([v, t]) => `<button role="tab" aria-selected="${this.view === v}" class="${this.view === v ? 'on' : ''}" data-c="view" data-v="${v}">${t}</button>`).join('')}</div>
<div class="tp-body">${this.view === 'tamano' ? this.sizeHTML(rbs) : this.view === 'gravedad' ? this.gravHTML(rbs) : this.tableHTML(rbs)}</div>`}`;
this.renderHits();
},
renderHits() {
const q = $('#cmp-q'), ul = $('#cmp-hits'); if (!q || !ul) return;
const s = norm(q.value.trim()); if (!s) { ul.innerHTML = ''; return; }
const hits = World.rb.filter(rb => this.can(rb) && !this.has(rb.id) && [rb.def.name, rb.def.short || '', ...(rb.def.aka || [])].some(n => norm(n).includes(s))).slice(0, 7);
ul.innerHTML = hits.map(rb => `<li><button data-c="add" data-id="${rb.id}"><i class="dot" style="--c:${rb.def.color}"></i>${esc(rb.def.name)}<small>${esc(TYPE_LABEL[rb.def.type])}</small></button></li>`).join('') || '<li class="note">Sin resultados.</li>';
},
tableHTML(rbs) {
const rows = this.ROWS.map(([label, re, k]) => {
const cells = rbs.map(rb => {
if (!re) return { txt: TYPE_LABEL[rb.def.type] + (rb.def.type === 'moon' && rb.parent ? ' de ' + rb.parent.def.name : '') };
if (k === 'au' && rb.def.type === 'moon') return { txt: 'No aplica: orbita a ' + rb.parent.def.name, na: true };
const f = this.field(rb, re); return f ? { txt: String(f.val), calc: f.kind === 'calc', n: k ? this.num(rb, k) : null } : { txt: 'No disponible', na: true };
});
if (re && cells.every(c => c.na)) return '';                        // campo sin sentido para estos objetos
const B = k ? this.bars(cells.map(c => c.n)) : null;
return `<tr><th scope="row">${label}${B && B.log ? '<small class="cmp-log" data-tip="Las diferencias son tan grandes que las barras usan escala logarítmica">escala log.</small>' : ''}</th>${cells.map((c, i) => `<td class="${c.na ? 'na' : ''}">${esc(c.txt)}${c.calc ? ' <small class="calc" data-tip="Valor calculado a partir de otros datos">calc.</small>' : ''}${B ? `<span class="cmp-bar"><i style="--c:${rbs[i].def.color};transform:scaleX(${B.k(c.n)})"></i></span>` : ''}</td>`).join('')}</tr>`;
}).join('');
return `<div class="cmp-scroll"><table class="cmp-table"><thead><tr><th></th>${rbs.map(rb => `<th scope="col">${esc(rb.def.short || rb.def.name)}</th>`).join('')}</tr></thead><tbody>${rows}</tbody></table></div>
<p class="note">Datos de referencia de las fichas de SOLARIS (NASA y JPL). «calc.» indica valores calculados a partir de la masa y el radio.${rbs.some(rb => rb.def.type === 'moon') ? ' En las lunas, el periodo y la velocidad orbital se refieren a su órbita alrededor del planeta.' : ''}</p>`;
},
sizeHTML(rbs) {
const D = rbs.map(rb => 2 * rb.def.R), W = 760, H = 310, pad = 18, gap = 18, SLOT = 112;   // espacio mínimo por etiqueta
const sum = D.reduce((a, c) => a + c, 0), s = Math.min(Math.max(1e-9, W - pad * 2 - gap * (rbs.length - 1) - SLOT * rbs.length * 0.5) / sum, (H - 84) / Math.max(...D));
let x = pad; const base = H - 58;
const circ = rbs.map((rb, i) => {
const r = D[i] * s / 2, w = Math.max(2 * r, SLOT), cx = x + w / 2, tiny = r < 1.5; x += w + gap;
const id = 'cg' + i;
return `<defs><radialGradient id="${id}" cx="35%" cy="35%" r="70%"><stop offset="0" stop-color="${rb.def.color}"/><stop offset="1" stop-color="${rb.def.color}" stop-opacity=".35"/></radialGradient></defs>
<circle cx="${cx}" cy="${base - Math.max(r, 1.5)}" r="${Math.max(r, 1.5)}" fill="url(#${id})"${tiny ? ' stroke="#f2c879" stroke-width="1"' : ''}/>
<text x="${cx}" y="${base + 18}" text-anchor="middle" class="sz-name">${esc(rb.def.short || rb.def.name)}</text>
<text x="${cx}" y="${base + 33}" text-anchor="middle" class="sz-val">${fmt(D[i], D[i] < 100 ? 1 : 0)} km</text>${tiny ? `<text x="${cx}" y="${base + 47}" text-anchor="middle" class="sz-val">casi invisible a esta escala</text>` : ''}`;
}).join('');
const vbW = Math.max(W, x - gap + pad);
const big = rbs[D.indexOf(Math.max(...D))], small = rbs[D.indexOf(Math.min(...D))], ratio = Math.max(...D) / Math.min(...D);
return `<svg class="cmp-size" viewBox="0 0 ${vbW} ${H}" role="img" aria-label="Comparación de tamaños a escala relativa">${circ}<line x1="${pad}" y1="${base}" x2="${vbW - pad}" y2="${base}" class="sz-base"/></svg>
<p class="cmp-ratio">El diámetro de <b>${esc(big.def.name)}</b> es <b>${fmt(ratio, ratio < 10 ? 1 : 0)} veces</b> el de <b>${esc(small.def.name)}</b>.</p>
<p class="note">Escala relativa basada en diámetro real. No usa la escala visual ampliada del Sistema Solar.</p>`;
},
gravHTML(rbs) {
const G = rbs.map(rb => this.grav(rb)), B = this.bars(G);
return `<div class="cmp-grav">${rbs.map((rb, i) => {
const g = G[i];
if (g == null) return `<div class="cg-row"><b>${esc(rb.def.name)}</b><span class="na">No disponible</span></div>`;
const f = g / 9.81, kg = 70 * f;
return `<div class="cg-row"><b><i class="dot" style="--c:${rb.def.color}"></i>${esc(rb.def.name)}</b><span class="cg-g">${fmt(g, g < 0.01 ? 5 : g < 1 ? 3 : 2)} m/s²</span>
<span class="cg-f">${fmt(f, f < 0.01 ? 4 : f < 1 ? 2 : 2)} × gravedad terrestre</span>${B ? `<span class="cmp-bar"><i style="--c:${rb.def.color};transform:scaleX(${B.k(g)})"></i></span>` : ''}
<small>Una báscula que marca 70 kg en la Tierra marcaría ≈ ${fmt(kg, kg < 1 ? 2 : 1)} kg${rb.def.type === 'planet' && ['jupiter', 'saturno', 'urano', 'neptuno'].includes(rb.id) ? ' (en la parte alta de la atmósfera: no tiene superficie sólida)' : ''}.</small></div>`;
}).join('')}</div>${B && B.log ? '<p class="note">Las barras usan escala logarítmica por la gran diferencia entre valores.</p>' : ''}
<p class="note">Factor respecto a la gravedad terrestre (9.81 m/s²). La masa no cambia: lo que varía es el peso.</p>`;
},
};
document.addEventListener('click', e => {
if (e.target.closest('#info-cmp')) { const id = e.target.closest('#info-cmp').dataset.cmp; Compare.has(id) ? Compare.remove(id) : Compare.add(id); }
else if (e.target.closest('[data-cmp-open]')) Compare.open();
else if (e.target.closest('[data-cmp-clear]')) Compare.clear();
});
const TL = {
MIN: Date.UTC(1950, 0, 1), MAX: Date.UTC(2130, 11, 31), cats: { exploracion: true, astronomia: true }, sel: null,
isOpen() { const p = $('#tl'); return !!p && p.classList.contains('open'); },
toggle(v) {
let p = $('#tl');
if (!p) { p = this.build(); }
const open = v == null ? !this.isOpen() : v;
p.classList.toggle('open', open); $('#app').classList.toggle('tl-open', open);
const b = $('#t-tl'); if (b) b.classList.toggle('on', open);
if (open) { UI.openSection(null); this.refresh(true); }
},
build() {
const p = document.createElement('section'); p.id = 'tl'; p.className = 'ui'; p.setAttribute('aria-label', 'Línea del tiempo');
const J = [['-10 años', -10, 'y'], ['-1 año', -1, 'y'], ['-1 mes', -1, 'm'], ['-1 día', -1, 'd'], ['+1 día', 1, 'd'], ['+1 mes', 1, 'm'], ['+1 año', 1, 'y'], ['+10 años', 10, 'y']];
p.innerHTML = `
<div class="tl-top">
<div class="tl-ind"><span class="tl-tag" id="tl-tag">SIMULACIÓN</span><b id="tl-date"></b><small id="tl-rel"></small></div>
<div class="tl-jumps" role="group" aria-label="Saltos de fecha">${J.map(([t, n, u]) => `<button class="tl-j" data-j="${n}" data-u="${u}">${t}</button>`).join('')}</div>
<div class="tl-right"><input id="tl-pick" type="date" min="1600-01-01" max="2400-12-31" aria-label="Elegir fecha"><button class="txt-btn" id="tl-today">Volver a hoy</button><button class="icon-btn" id="tl-x" aria-label="Cerrar la línea del tiempo">&times;</button></div>
</div>
<div class="tl-track-wrap"><div class="tl-ticks" id="tl-ticks" aria-hidden="true"></div><div class="tl-evs" id="tl-evs"></div>
<input id="tl-range" type="range" min="${this.MIN}" max="${this.MAX}" step="86400000" aria-label="Desplazar la fecha de la simulación"></div>
<div class="tl-bottom">
<div class="tl-cats" role="group" aria-label="Categorías de eventos">${Object.entries(EVENT_CATS).map(([k, t]) => `<label class="tl-cat tl-${k}"><input type="checkbox" data-cat="${k}" checked><i></i>${t}</label>`).join('')}</div>
<div class="tl-card" id="tl-card" hidden></div>
</div>`;
$('#app').appendChild(p);
const span = this.MAX - this.MIN, years = [];
for (let y = 1950; y <= 2130; y += 10) years.push(`<span style="left:${(Date.UTC(y, 0, 1) - this.MIN) / span * 100}%">${y}</span>`);
$('#tl-ticks').innerHTML = years.join('');
this.paintEvents();
p.addEventListener('click', e => {
const j = e.target.closest('[data-j]');
if (j) { const d = jdToDate(Time.jd), n = +j.dataset.j, u = j.dataset.u;
if (u === 'd') d.setUTCDate(d.getUTCDate() + n); else if (u === 'm') d.setUTCMonth(d.getUTCMonth() + n); else d.setUTCFullYear(d.getUTCFullYear() + n);
UI.setSimDate(dateToJD(d), { toast: false }); }
const ev = e.target.closest('[data-ev]'); if (ev) this.showEvent(+ev.dataset.ev);
const a = e.target.closest('[data-tla]');
if (a) { const E = EVENTS[this.sel], rb = World.byId[E.obj];
if (a.dataset.tla === 'date') { UI.setSimDate(dateToJD(new Date(E.date + 'T12:00:00Z')), { pause: true }); if (rb && !rb.hidden) UI.select(rb, { fly: true }); }
else if (a.dataset.tla === 'go' && rb) { if (rb.hidden) UI.toast(rb.def.name + ' no existe en la fecha simulada.'); else UI.select(rb, { fly: true }); } }
});
$('#tl-x').addEventListener('click', () => this.toggle(false));
$('#tl-today').addEventListener('click', () => { Time.goLive(); UI.refreshTime(); this.refresh(); UI.toast('De vuelta a la fecha y hora actuales'); });
$('#tl-pick').addEventListener('change', e => { if (e.target.value) UI.setSimDate(Astro.jdFromISO(e.target.value) + 0.5); });
p.addEventListener('change', e => { const c = e.target.closest('[data-cat]'); if (c) { this.cats[c.dataset.cat] = c.checked; this.paintEvents(); } });
const r = $('#tl-range');
r.addEventListener('input', () => { this.pending = +r.value; Time.paused = true; if (!this.raf) this.raf = requestAnimationFrame(() => { this.raf = 0; UI.setSimDate(dateToJD(new Date(this.pending)), { toast: false }); }); });
r.addEventListener('change', () => { const y = new Date(+r.value).getUTCFullYear(); if (y < 1800 || y > 2050) UI.toast('Fuera de 1800–2050 la precisión de las posiciones disminuye.'); });
return p;
},
visibleEvents() { return EVENTS.map((e, i) => ({ e, i })).filter(({ e }) => this.cats[e.cat]); },
paintEvents() {
const span = this.MAX - this.MIN;
$('#tl-evs').innerHTML = this.visibleEvents().map(({ e, i }) => { const t = Date.parse(e.date + 'T12:00:00Z'); if (t < this.MIN || t > this.MAX) return '';
return `<button class="tl-ev tl-${e.cat}${e.pred ? ' pred' : ''}${this.sel === i ? ' on' : ''}" data-ev="${i}" style="left:${(t - this.MIN) / span * 100}%" data-tip="${esc(e.title)} · ${this.fmtDate(e.date)}" aria-label="${esc(e.title)}, ${this.fmtDate(e.date)}"></button>`; }).join('');
},
fmtDate(iso) { return new Date(iso + 'T12:00:00Z').toLocaleDateString(I18N.loc(), { day: 'numeric', month: 'long', year: 'numeric', timeZone: 'UTC' }); },
showEvent(i) {
this.sel = i; const E = EVENTS[i], rb = World.byId[E.obj];
const c = $('#tl-card'); c.hidden = false;
c.innerHTML = `<span class="tl-cat-tag tl-${E.cat}">${EVENT_CATS[E.cat]}${E.pred ? ' · fecha prevista' : ''}</span><h4>${esc(E.title)}</h4><time>${this.fmtDate(E.date)}</time><p>${esc(E.desc)}</p>
<dl>${rb ? `<dt>Objeto relacionado</dt><dd>${esc(rb.def.name)}</dd>` : ''}${E.mission ? `<dt>Misión</dt><dd>${esc(E.mission)}</dd>` : ''}</dl>
<div class="tl-card-act"><button class="cta" data-tla="date">Ver fecha</button>${rb ? '<button class="cta-ghost" data-tla="go">Ir al objeto</button>' : ''}</div>`;
this.paintEvents();
},
refresh(force) {
const sim = UI.isSimulatedDate(), d = jdToDate(Time.jd);
const chip = $('#t-sim'); if (chip) chip.hidden = !sim;
if (!this.isOpen() && !force) return;
$('#tl-date').textContent = d.toLocaleDateString(I18N.loc(), { day: '2-digit', month: 'short', year: 'numeric', timeZone: 'UTC' }).replace('.', '').toUpperCase();
$('#tl-tag').textContent = sim ? 'SIMULACIÓN' : 'EN VIVO'; $('#tl-tag').classList.toggle('live', !sim);
const days = Time.jd - dateToJD(new Date()), ad = Math.abs(days);
$('#tl-rel').textContent = !sim ? 'Fecha y hora actuales' : (ad < 60 ? fmt(ad, 0) + ' días' : ad < 730 ? fmt(ad / 30.44, 0) + ' meses' : fmt(ad / 365.25, 0) + ' años') + (days > 0 ? ' en el futuro' : ' en el pasado');
const r = $('#tl-range'); if (document.activeElement !== r) r.value = clamp(d.getTime(), this.MIN, this.MAX);
const pk = $('#tl-pick'); if (document.activeElement !== pk) pk.value = d.toISOString().slice(0, 10);
},
};
$('#t-tl').addEventListener('click', () => TL.toggle());
const SCI_LAYERS = [
{ k: 'hz', name: 'Zona habitable', desc: 'Región aproximada donde podría existir agua líquida en la superficie de un planeta' },
];
const SCI_PRESETS = {
exploracion: { name: 'Exploración', desc: 'Interfaz limpia con pocas capas', L: { orbits: true, moonOrbits: true, namesPlanets: true, namesMoons: false, namesSmall: false, orbitInfo: false, cometOrbits: false, hz: false } },
educativo: { name: 'Educativo', desc: 'Órbitas, nombres y zona habitable', L: { orbits: true, moonOrbits: true, namesPlanets: true, namesMoons: true, namesSmall: false, orbitInfo: false, hz: true } },
astronomia: { name: 'Astronomía', desc: 'Órbitas, trayectorias, datos orbitales y zona habitable', L: { orbits: true, moonOrbits: true, cometOrbits: true, namesPlanets: true, namesMoons: true, orbitInfo: true, hz: true } },
cinematico: { name: 'Cinemático', desc: 'Sin capas técnicas', L: { orbits: false, moonOrbits: false, cometOrbits: false, namesPlanets: false, namesMoons: false, namesSmall: false, namesCraft: false, orbitInfo: false, markers: false, constLines: false, constNames: false, hz: false } },
};
const AU_MKM = AU_KM / 1e6;
const Sci = {
hz() { return HZ_MODELS[S.hzModel] || HZ_MODELS.conservador; },
status(r) {
const H = this.hz();
if (r < H.in * 0.9) return ['fuera', 'Fuera (más cerca del Sol)'];
if (r < H.in) return ['borde', 'Cerca del límite interior, por fuera'];
if (r < H.in * 1.1) return ['dentro', 'Dentro, cerca del límite interior'];
if (r <= H.out * 0.9) return ['dentro', 'Dentro'];
if (r <= H.out) return ['dentro', 'Dentro, cerca del límite exterior'];
if (r <= H.out * 1.1) return ['borde', 'Cerca del límite exterior, por fuera'];
return ['fuera', 'Fuera (más lejos del Sol)'];
},
setModel(k) { S.hzModel = k; Settings.state.simulation.hzModel = k; Settings.save(); this.paintPop(); this.paintSection(); },
applyPreset(id) {
const P = SCI_PRESETS[id]; if (!P) return;
for (const k in P.L) if (S.layers[k] !== P.L[k]) UI.setLayer(k, P.L[k]);
Settings.state.simulation.sciPreset = id; Settings.save();
this.paintSection(); UI.toast('Preset de capas: ' + P.name);
},
sectionHTML() {
const pr = Settings.state.simulation.sciPreset, H = this.hz();
return `<h3 class="sci-h">Capas científicas</h3>
<div class="seg sci-presets" role="radiogroup" aria-label="Presets de capas">${Object.entries(SCI_PRESETS).map(([k, p]) => `<button role="radio" aria-checked="${pr === k}" class="${pr === k ? 'on' : ''}" data-sci-preset="${k}" data-tip="${p.desc}">${p.name}</button>`).join('')}</div>
${SCI_LAYERS.map(l => `<div class="sci-row"><label class="mix-sw sci-sw"><input type="checkbox" data-layer="${l.k}" ${S.layers[l.k] ? 'checked' : ''}><span class="sw" aria-hidden="true"></span><span>${l.name}</span></label>
<button class="sci-q" data-sci-info="${l.k}" aria-label="Información sobre: ${l.name}" data-tip="Información">?</button></div>
<p class="sci-state" data-sci-state="${l.k}">${S.layers[l.k] ? 'Activada' : 'Desactivada'} · modelo ${H.name.toLowerCase()} (${fmt(H.in, 2)}–${fmt(H.out, 2)} UA)</p>`).join('')}
<p class="note">Visualizaciones educativas superpuestas: no modifican las órbitas ni los datos de los objetos. En Modo Cine se ocultan y se restauran al salir.</p>`;
},
paintSection() {
const st = document.querySelector('[data-sci-state="hz"]'), H = this.hz();
if (st) st.textContent = `${S.layers.hz ? 'Activada' : 'Desactivada'} · modelo ${H.name.toLowerCase()} (${fmt(H.in, 2)}–${fmt(H.out, 2)} UA)`;
const pr = Settings.state.simulation.sciPreset;
document.querySelectorAll('[data-sci-preset]').forEach(b => { const on = b.dataset.sciPreset === pr; b.classList.toggle('on', on); b.setAttribute('aria-checked', on); });
},
openPop() {
let p = $('#sci-pop');
if (!p) {
p = document.createElement('section'); p.id = 'sci-pop'; p.className = 'ui'; p.setAttribute('role', 'dialog'); p.setAttribute('aria-labelledby', 'sci-pop-t');
$('#app').appendChild(p);
p.addEventListener('click', e => { if (e.target.closest('[data-sci-close]')) this.closePop(); const m = e.target.closest('[data-hz-model]'); if (m) this.setModel(m.dataset.hzModel); });
addEventListener('keydown', e => { if (e.key === 'Escape' && p.classList.contains('open')) { e.stopImmediatePropagation(); this.closePop(); } }, true);
}
p.classList.add('open'); this.paintPop();
clearInterval(this.popT); this.popT = setInterval(() => this.paintPlanets(), 500);   // las posiciones cambian con la fecha
},
closePop() { const p = $('#sci-pop'); if (p) p.classList.remove('open'); clearInterval(this.popT); },
paintPop() {
const p = $('#sci-pop'); if (!p || !p.classList.contains('open')) return;
const H = this.hz(), lim = (d, n) => `<b>${fmt(d, 2)} UA</b> · ${fmt(d * AU_MKM, 1)} millones de km <small>(${n})</small>`;
const scaleTxt = S.scale === 'real' ? 'Escala real: el anillo y las órbitas usan distancias reales proporcionales.'
: 'Escala ' + String((SCALES[S.scale] && SCALES[S.scale].name) || 'visual').toLowerCase() + ': las distancias se comprimen para facilitar la exploración. El anillo usa exactamente la misma transformación que las órbitas, así que la posición de los planetas respecto a la zona es correcta, aunque las proporciones entre distancias no sean las reales.';
p.innerHTML = `
<header class="tp-head"><div><span class="tp-kicker">Capa científica</span><h2 id="sci-pop-t">Zona habitable</h2></div><button class="icon-btn" data-sci-close aria-label="Cerrar">&times;</button></header>
<p>Región aproximada alrededor de una estrella donde podrían existir temperaturas compatibles con agua líquida en la superficie de un planeta, dependiendo de factores como su atmósfera, presión y composición.</p>
<p class="sci-warn">Estar dentro de la zona habitable no significa automáticamente que un planeta sea habitable.</p>
<h4>Modelo</h4>
<div class="seg" role="radiogroup" aria-label="Modelo de zona habitable">${Object.entries(HZ_MODELS).map(([k, m]) => `<button role="radio" aria-checked="${S.hzModel === k}" class="${S.hzModel === k ? 'on' : ''}" data-hz-model="${k}">${m.name}</button>`).join('')}</div>
<dl class="sci-lims"><dt>Límite interior</dt><dd>${lim(H.in, H.inName)}</dd><dt>Límite exterior</dt><dd>${lim(H.out, H.outName)}</dd></dl>
<p class="note">Fuente: ${H.ref}.${H.calc ? ' Valores calculados a partir de los datos publicados.' : ''} Los límites son aproximados y dependen del modelo climático utilizado.</p>
<h4>Planetas respecto a la zona <small id="sci-date"></small></h4>
<div id="sci-planets" class="sci-planets"></div>
<p class="note">${scaleTxt}</p>`;
this.paintPlanets();
},
paintPlanets() {
const el = $('#sci-planets'); if (!el) return;
const pl = World.rb.filter(r => r.def.type === 'planet');
el.innerHTML = pl.map(rb => { const r = V.len(rb.helio), [c, t] = this.status(r);
return `<div class="sp-row sp-${c}${rb.id === 'tierra' ? ' sp-ref' : ''}"><i class="dot" style="--c:${rb.def.color}"></i><b>${esc(rb.def.name)}</b><span>${fmt(r, 2)} UA</span><em>${t}</em></div>`; }).join('');
const d = $('#sci-date'); if (d) d.textContent = '· ' + UI.dateStr(Time.jd);
},
labels() {
let A = $('#hz-labels');
if (!A) { A = document.createElement('div'); A.id = 'hz-labels'; A.innerHTML = '<span class="hz-lab" id="hz-in"></span><span class="hz-lab" id="hz-out"></span>'; $('#labels').appendChild(A); }
const show = S.layers.hz && (World.cineK || 0) < 0.4 && !Flight.on;
A.hidden = !show; if (!show) return;
const H = this.hz(), cam = World.cam, h = Math.hypot(cam[0], cam[2]) || 1, dir = [cam[0] / h, 0, cam[2] / h];
const put = (id, au, name) => {
const el = $('#' + id), r = ScaleState.dist(au), o = World.project([dir[0] * r, 0, dir[2] * r], {}), o2 = World.project([-dir[0] * r, 0, -dir[2] * r], {});
const big = o.on && o2.on && Math.hypot(o.x - o2.x, o.y - o2.y) > 90;           // solo si el anillo se ve con tamaño suficiente
el.hidden = !(o.on && big); if (el.hidden) return;
el.style.transform = `translate(${o.x}px, ${o.y}px)`;
const t = `${name} · ${fmt(au, 2)} UA`; if (el.textContent !== t) { el.textContent = t; el.dataset.tip = `${name}: ${fmt(au, 2)} UA ≈ ${fmt(au * AU_MKM, 1)} millones de km del Sol (modelo ${H.name.toLowerCase()})`; }
};
put('hz-in', H.in, 'Límite interior'); put('hz-out', H.out, 'Límite exterior');
},
};
const _secCapas = UI.sectionHTML;
UI.sectionHTML = function (id) { const h = _secCapas.call(this, id); return id === 'capas' ? h + Sci.sectionHTML() : h; };
UI.TOOL_ROWS = (UI.TOOL_ROWS || []).concat([['science', ICON.capas || ICON.herramientas, 'Capas científicas', 'Zona habitable y presets de capas']]);
UI.TOOL_ACTIONS = Object.assign(UI.TOOL_ACTIONS || {}, { science: () => UI.openSection('capas') });
document.addEventListener('click', e => {
const q = e.target.closest('[data-sci-info]'); if (q) Sci.openPop();
const pr = e.target.closest('[data-sci-preset]'); if (pr) Sci.applyPreset(pr.dataset.sciPreset);
});
const _setLayer = UI.setLayer;
UI.setLayer = function (k, v) {
_setLayer.call(this, k, v);
if (k === 'hz') Sci.paintSection();
if (!Sci.applying && Settings.state.simulation.sciPreset) { const P = SCI_PRESETS[Settings.state.simulation.sciPreset]; if (P && k in P.L && P.L[k] !== v) { Settings.state.simulation.sciPreset = null; Settings.save(); Sci.paintSection(); } }
};
const _apply = Sci.applyPreset.bind(Sci); Sci.applyPreset = id => { Sci.applying = true; try { _apply(id); } finally { Sci.applying = false; } };
const MU_SUN = G_CONST * BODY.sol.mass;              // m³/s², a partir de la masa del Sol de los datos
const AU_M = AU_KM * 1000;
const PLAN_GROUPS = [['Sol y planetas', r => r.isSun || r.def.type === 'planet'], ['Lunas', r => r.def.type === 'moon'], ['Planetas enanos y transneptunianos', r => r.def.type === 'dwarf' || r.def.type === 'tno'],
['Asteroides', r => r.def.type === 'asteroid'], ['Cometas', r => r.def.type === 'comet'], ['Naves y sondas', r => r.isCraft]];
const Planner = {
origin: 'tierra', dest: 'marte', ship: null, route: null, routeOn: false,
get mode() { return Settings.state.simulation.plannerMode || 'simple'; },
set mode(v) { Settings.state.simulation.plannerMode = v; Settings.save(); },
isOpen() { const p = $('#plan'); return !!p && p.classList.contains('open'); },
open(opts) {
opts = opts || {};
if (opts.origin) this.origin = opts.origin; if (opts.dest) this.dest = opts.dest;
let p = $('#plan');
if (!p) {
p = document.createElement('section'); p.id = 'plan'; p.className = 'tool-panel ui'; p.setAttribute('role', 'dialog'); p.setAttribute('aria-labelledby', 'plan-title'); $('#app').appendChild(p);
p.addEventListener('change', e => {
const t = e.target;
if (t.id === 'pl-o') this.origin = t.value; else if (t.id === 'pl-d') this.dest = t.value; else if (t.id === 'pl-s') this.ship = t.value; else return;
if (this.routeOn) this.buildRoute(); this.render();
});
p.addEventListener('click', e => {
const a = e.target.closest('[data-pl]'); if (!a) return; const v = a.dataset.pl;
if (v === 'close') this.close();
else if (v === 'mode') { this.mode = a.dataset.m; if (this.routeOn) this.buildRoute(); this.render(); }
else if (v === 'swap') { [this.origin, this.dest] = [this.dest, this.origin]; if (this.routeOn) this.buildRoute(); this.render(); }
else if (v === 'see-o' || v === 'see-d') { this.close(); UI.select(World.byId[v === 'see-o' ? this.origin : this.dest], { fly: true }); }
else if (v === 'route') { this.routeOn ? this.hideRoute() : this.showRoute(); this.render(); }
else if (v === 'window') { const H = this.calc().hoh; if (H && H.win != null) { UI.setSimDate(Time.jd + H.win, { pause: true }); World.update(Time.jd, 0, 0); if (this.routeOn) this.buildRoute(); this.render(); } }   // posiciones de la nueva fecha antes de recalcular
else if (v === 'start') this.start();
});
addEventListener('keydown', e => { if (e.key === 'Escape' && this.isOpen()) { e.stopImmediatePropagation(); this.close(); } }, true);
}
UI.openSection(null); Compare.close && Compare.close();
p.classList.add('open'); this.render();
clearInterval(this.tick); this.tick = setInterval(() => { if (this.isOpen()) this.paintResults(); }, 1000);   // la fecha global puede cambiar
},
close() { const p = $('#plan'); if (p) p.classList.remove('open'); clearInterval(this.tick); },
rb(id) { const r = World.byId[id]; return r && !r.hidden ? r : null; },
shipDef() { return SHIPS.find(s => s.id === this.ship) || SHIPS.find(s => s.id === 'ranger') || SHIPS[0]; },
orbitalOK(o, d) { return o && d && o !== d && !o.isSun && !d.isSun && !o.isCraft && !d.isCraft && o.parent && o.parent.isSun && d.parent && d.parent.isSun; },
meanAU(rb) { const n = Compare.num(rb, 'au'); return n || V.len(rb.helio); },
calc() {
const o = this.rb(this.origin), d = this.rb(this.dest), sh = this.shipDef(); if (!o || !d || o === d) return { o, d, sh };
const dAU = V.len(V.sub(d.helio, o.helio)), km = dAU * AU_KM;
const days = Flight.tripDays(sh, km);
const R = { o, d, sh, dAU, km, days, arr: Time.jd + days };
if (this.mode === 'orbital' && this.orbitalOK(o, d)) {
const r1 = this.meanAU(o), r2 = this.meanAU(d), a = (r1 + r2) / 2;
const T = Math.PI * Math.sqrt(Math.pow(a * AU_M, 3) / MU_SUN) / 86400;                   // días
const v = r => Math.sqrt(MU_SUN / (r * AU_M)), vt = r => Math.sqrt(MU_SUN * (2 / (r * AU_M) - 1 / (a * AU_M)));
const dv1 = Math.abs(vt(r1) - v(r1)) / 1000, dv2 = Math.abs(v(r2) - vt(r2)) / 1000;
const P1 = Info.periodDays(o), P2 = Info.periodDays(d);
let win = null, phase = null, phi = null, syn = null;
if (P1 && P2) {
const n1 = 2 * Math.PI / P1.v, n2 = 2 * Math.PI / P2.v, w = n2 - n1, lon = rb => Math.atan2(rb.helio[1], rb.helio[0]);
const TAU = Math.PI * 2, wrap = x => ((x % TAU) + TAU) % TAU;
phase = wrap(lon(d) - lon(o)); phi = wrap(Math.PI - n2 * T);
win = Math.abs(w) < 1e-9 ? null : w > 0 ? wrap(phi - phase) / w : wrap(phase - phi) / -w;
const dphi = Math.abs(((phase - phi + Math.PI * 3) % TAU) - Math.PI);            // margen propio de usar movimientos medios
if (dphi < 6 * Math.PI / 180) win = 0;
syn = Math.abs(w) < 1e-9 ? null : TAU / Math.abs(w);
}
R.hoh = { r1, r2, a, T, dv1, dv2, win, phase, phi, syn };
const vel = rb => { const o2 = rb.def.orbit; if (!o2 || !Astro.helioPos) return null; try { const p0 = Astro.helioPos(o2, Time.jd - 0.5), p1 = Astro.helioPos(o2, Time.jd + 0.5); return V.sub(p1, p0); } catch (e) { return null; } };
const vo = vel(o), vd = vel(d); R.vrel = vo && vd ? V.len(V.sub(vd, vo)) * AU_KM / 86400 : null;
}
return R;
},
optionsHTML(sel) {
return PLAN_GROUPS.map(([t, f]) => { const it = World.rb.filter(r => f(r) && !r.hidden); return it.length ? `<optgroup label="${t}">${it.map(r => `<option value="${r.id}" ${r.id === sel ? 'selected' : ''}>${esc(r.def.name)}</option>`).join('')}</optgroup>` : ''; }).join('');
},
render() {
const p = $('#plan'); if (!p) return;
if (!this.ship) this.ship = this.shipDef().id;
this._lastHTML = null;
p.innerHTML = `
<header class="tp-head"><div><span class="tp-kicker">Herramienta</span><h2 id="plan-title">Planificador de viaje</h2></div><button class="icon-btn" data-pl="close" aria-label="Cerrar">&times;</button></header>
<div class="seg tp-tabs" role="radiogroup" aria-label="Modo de planificación">${[['simple', 'Modo simple'], ['orbital', 'Modo orbital']].map(([m, t]) => `<button role="radio" aria-checked="${this.mode === m}" class="${this.mode === m ? 'on' : ''}" data-pl="mode" data-m="${m}">${t}</button>`).join('')}</div>
<div class="pl-form">
<label>Origen<select id="pl-o">${this.optionsHTML(this.origin)}</select></label>
<button class="icon-btn pl-swap" data-pl="swap" aria-label="Intercambiar origen y destino" data-tip="Intercambiar">⇄</button>
<label>Destino<select id="pl-d">${this.optionsHTML(this.dest)}</select></label>
<label>Nave<select id="pl-s">${SHIPS.map(s => `<option value="${s.id}" ${s.id === this.ship ? 'selected' : ''}>${esc(s.name)}</option>`).join('')}</select></label>
</div>
<div id="pl-res"></div>
<button class="txt-btn pl-win" id="pl-winb" data-pl="window" hidden>Ir a la fecha de la ventana de transferencia</button>
<div class="pl-act">
<button class="cta-ghost" data-pl="see-o">Ver origen</button><button class="cta-ghost" data-pl="see-d">Ver destino</button>
<button class="cta-ghost${this.routeOn ? ' on' : ''}" data-pl="route">${this.routeOn ? 'Ocultar trayectoria' : 'Visualizar viaje'}</button>
<button class="cta" data-pl="start">Iniciar viaje</button>
</div>`;
this.paintResults();
},
paintResults() {
const el = $('#pl-res'); if (!el) return;
const R = this.calc(), sh = R.sh;
if (!R.o || !R.d) { this._lastHTML = null; el.innerHTML = '<p class="note">El origen o el destino no existe en la fecha simulada (por ejemplo, una nave aún no lanzada). Cambia la fecha o elige otro objeto.</p>'; return; }
if (R.o === R.d) { this._lastHTML = null; el.innerHTML = '<p class="note">Elige un destino distinto del origen.</p>'; return; }
const date = jd => UI.dateStr(jd), dist = km => km >= 1e6 ? fmt(km / 1e6, km >= 1e8 ? 0 : 1) + ' millones de km' : fmt(km, 0) + ' km';
const fict = sh.info && sh.info.franquicia ? `Tecnología ficticia (${esc(sh.info.franquicia)})` : 'Nave conceptual de SOLARIS';
let html = `<div class="pl-grid">
<div><span>Distancia actual</span><b>${dist(R.km)}</b><small>${fmt(R.km, 0)} km · ${fmt(R.dAU, R.dAU < 0.01 ? 5 : 3)} UA</small></div>
<div><span>Nave</span><b>${esc(sh.name)}</b><small>Velocidad máxima ${fmt(sh.maxSpeed, 0)} km/s · aceleración ${fmt(sh.accel / 9.81, 1)} g</small></div>
<div><span>Tiempo estimado</span><b>${fmtDur(R.days)}</b><small>Acelerando, en crucero y frenando en línea recta</small></div>
<div><span>Salida</span><b>${date(Time.jd)}</b><small>Fecha de la simulación</small></div>
<div><span>Llegada aproximada</span><b>${date(R.arr)}</b><small>${UI.isSimulatedDate() ? 'Según la fecha simulada' : 'Saliendo ahora'}</small></div>
</div>
<p class="pl-note">${fict} — cálculo realizado con los parámetros internos de SOLARIS. Estimación simplificada basada en distancia y velocidad seleccionada. No representa una trayectoria orbital real${R.d.parent && !R.d.isSun ? ' ni el movimiento del destino durante el viaje' : ''}.</p>`;
if (this.mode === 'orbital') {
if (!R.hoh) html += `<p class="pl-note">El modo orbital necesita dos cuerpos que orbiten directamente al Sol (por ejemplo, Tierra → Marte). Para lunas, naves o el propio Sol se muestra solo la estimación simple.</p>`;
else {
const H = R.hoh, deg = x => fmt(x * 180 / Math.PI, 0) + '°';
html += `<h4 class="pl-h">Transferencia orbital simplificada (Hohmann)</h4><div class="pl-grid">
<div><span>Duración de la transferencia</span><b>${fmtDur(H.T)}</b><small>Media elipse entre ${fmt(H.r1, 2)} y ${fmt(H.r2, 2)} UA</small></div>
<div><span>Δv heliocéntrico</span><b>${fmt(H.dv1 + H.dv2, 1)} km/s</b><small>${fmt(H.dv1, 2)} al salir + ${fmt(H.dv2, 2)} al llegar</small></div>
<div><span>Próxima ventana de transferencia</span><b>${H.win != null ? date(Time.jd + H.win) : 'No disponible'}</b><small>${H.win != null ? (H.win < 1 ? 'Ventana abierta ahora (dentro del margen de esta aproximación)' : 'Dentro de ' + fmtDur(H.win)) + (H.syn ? ' · se repite cada ' + fmtDur(H.syn) : '') : ''}</small></div>
<div><span>Ángulo de fase</span><b>${H.phase != null ? deg(H.phase) : 'No disponible'}</b><small>${H.phi != null ? 'Necesario para partir: ' + deg(H.phi) : ''}</small></div>
<div><span>Velocidad relativa actual</span><b>${R.vrel != null ? fmt(R.vrel, 1) + ' km/s' : 'No disponible'}</b><small>Entre ambos cuerpos alrededor del Sol</small></div>
</div>
<p class="pl-note">Simplificaciones: órbitas circulares y coplanarias con el radio medio de cada cuerpo; solo gravedad del Sol; el Δv no incluye escapar del planeta de origen ni la captura en el destino. Al visualizarla, la trayectoria parte de la posición del origen en la próxima ventana. Es una aproximación educativa, no un cálculo de misión.</p>`;
}
}
if (html !== this._lastHTML) { el.innerHTML = html; this._lastHTML = html; }          // solo si cambió (no recrea nodos sin necesidad)
const wb = $('#pl-winb'), want = !!(R.hoh && R.hoh.win != null && R.hoh.win >= 1); if (wb && wb.hidden === want) wb.hidden = !want;
},
buildRoute() {
const R = this.calc(); if (!R.o || !R.d || R.o === R.d) { this.route = null; return; }
if (R.hoh) {
const H = R.hoh, e = Math.abs(H.r2 - H.r1) / (H.r2 + H.r1), out = H.r2 > H.r1;
const lon0 = Math.atan2(R.o.helio[1], R.o.helio[0]) + (H.win || 0) * 2 * Math.PI / Info.periodDays(R.o).v, N = 160, pts = [];
for (let i = 0; i < N; i++) { const nu = (out ? 0 : Math.PI) + Math.PI * i / (N - 1), r = H.a * (1 - e * e) / (1 + e * Math.cos(nu)), ang = lon0 + (out ? nu : nu - Math.PI); pts.push([r * Math.cos(ang), r * Math.sin(ang), 0]); }
this.route = { kind: 'hohmann', pts };
} else this.route = { kind: 'line', o: R.o.id, d: R.d.id };
},
showRoute() {
this.buildRoute(); if (!this.route) return;
this.routeOn = true; this.close();
const o = World.byId[this.origin], d = World.byId[this.dest];
if (this.route.kind === 'hohmann') { const r = ScaleState.dist(Math.max(this.calc().hoh.r2, this.calc().hoh.r1)); Cam.travelTo([0, 0, 0], r * 3.2, 1.0); }
else { const c = V.scale(V.add(o.posS, d.posS), 0.5), s = V.dist(o.posS, d.posS); Cam.travelTo(c, Math.max(s * 1.9, Math.max(o.rS, d.rS) * 8), 0.75); }
UI.toast(this.route.kind === 'hohmann' ? 'Trayectoria de transferencia simplificada (Hohmann)' : 'Línea recta educativa: no representa la trayectoria física real');
},
hideRoute() { this.routeOn = false; this.route = null; const m = $('#pl-marks'); if (m) m.hidden = true; },
scenePts() {
const R = this.route; if (!R) return null;
if (R.kind === 'line') { const o = World.byId[R.o], d = World.byId[R.d], N = 64, out = [];
for (let i = 0; i < N; i++) out.push(V.add(o.posS, V.scale(V.sub(d.posS, o.posS), i / (N - 1)))); return out; }
return R.pts.map(p => ScaleState.mapVec(Astro.eclToScene(p)));
},
marks() {
let M = $('#pl-marks');
if (!M) { M = document.createElement('div'); M.id = 'pl-marks'; M.innerHTML = '<span class="pl-mk" id="pl-mo"></span><span class="pl-mk" id="pl-md"></span>'; $('#labels').appendChild(M); }
const show = this.routeOn && (World.cineK || 0) < 0.4 && !Flight.on; M.hidden = !show; if (!show) return;
[['pl-mo', this.origin, 'Origen'], ['pl-md', this.dest, 'Destino']].forEach(([id, bid, t]) => {
const rb = World.byId[bid], el = $('#' + id); if (!rb) { el.hidden = true; return; }
const o = World.project(rb.posS, {}); el.hidden = !o.on; if (!o.on) return;
el.style.transform = `translate(${o.x}px, ${o.y}px)`; const txt = t + ' · ' + rb.def.name; if (el.textContent !== txt) el.textContent = txt;
});
},
start() {
const R = this.calc(); if (!R.o || !R.d || R.o === R.d) return;
let sid = R.o.id;
if (!STARTS.some(s => s[0] === sid)) {       // el vuelo parte de una de las órbitas disponibles: la más cercana al origen
sid = STARTS.map(s => s[0]).filter(id => id !== 'luna' || (R.o.parent && R.o.parent.id === 'tierra')).sort((a, b) => V.len(V.sub(World.byId[a].helio, R.o.helio)) - V.len(V.sub(World.byId[b].helio, R.o.helio)))[0];
UI.toast('El vuelo comienza en la órbita disponible más cercana: ' + STARTS.find(s => s[0] === sid)[1]);
}
this.close(); this.hideRoute();
Flight.startId = sid; Flight.openHangar();
const i = SHIPS.findIndex(s => s.id === R.sh.id); if (i >= 0) Flight.setShip(i);
Flight.renderHangar(); Flight.launch();
Flight.setTarget(R.d); if (!Flight.intercept) Flight.toggle('intercept');
},
};
UI.TOOL_ROWS = (UI.TOOL_ROWS || []).concat([['planner', ICON.travel, 'Planificador de viaje', 'Distancia, tiempo y transferencia orbital entre dos objetos']]);
UI.TOOL_ACTIONS = Object.assign(UI.TOOL_ACTIONS || {}, { planner: () => Planner.open() });
document.addEventListener('click', e => { const b = e.target.closest('#info-plan'); if (b) Planner.open({ dest: b.dataset.plan }); });
ICON.gargantua = '<svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.5" aria-hidden="true"><circle cx="12" cy="12" r="3.2" fill="currentColor"/><ellipse cx="12" cy="12" rx="9.5" ry="3"/><path d="M5 9.5a7.5 7.5 0 0 1 14 0" opacity=".6"/></svg>';
ICON.systems = '<svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="1.5" aria-hidden="true"><circle cx="7" cy="12" r="2.6"/><circle cx="17" cy="12" r="2.6" fill="currentColor"/><path d="M9.8 12h4.4"/></svg>';
const SECTIONS_SOLAR = UI.SECTIONS.slice();
const SECTIONS_BY_SYSTEM = { gargantua: [['gargantua', 'Sistema Gargantúa'], ['capas', 'Capas'], ['config', 'Ajustes']] };
const Systems = {
busy: false,
current() { return World.system; },
openPicker() {
let p = $('#sys-pick');
if (!p) {
p = document.createElement('section'); p.id = 'sys-pick'; p.className = 'ui'; p.setAttribute('role', 'dialog'); p.setAttribute('aria-labelledby', 'sys-pick-t');
$('#app').appendChild(p);
p.addEventListener('click', e => { const c = e.target.closest('[data-sys]'); if (c) { this.closePicker(); this.go(c.dataset.sys); } if (e.target.closest('[data-sys-close]')) this.closePicker(); });
addEventListener('keydown', e => { if (e.key === 'Escape' && p.classList.contains('open')) { e.stopImmediatePropagation(); this.closePicker(); } }, true);
}
p.innerHTML = `<header class="tp-head"><div><span class="tp-kicker">Exploración</span><h2 id="sys-pick-t">Sistemas</h2></div><button class="icon-btn" data-sys-close aria-label="Cerrar">&times;</button></header>
<div class="sys-cards">${Object.values(SYSTEMS).map(s => `<button class="sys-card sys-${s.id}${s.id === World.system ? ' on' : ''}" data-sys="${s.id}" ${s.id === World.system ? 'aria-current="true"' : ''}>
<span class="sys-ico">${s.id === 'solar' ? ICON.overview : ICON.gargantua}</span><b>${s.name}</b><small>${s.sub}</small>${s.id === World.system ? '<em>Sistema actual</em>' : ''}</button>`).join('')}</div>`;
p.classList.add('open');
},
closePicker() { const p = $('#sys-pick'); if (p) p.classList.remove('open'); },
go(id) {
if (this.busy || id === World.system || !SYSTEMS[id]) return;
if (Flight.on) { UI.toast('Termina el vuelo antes de cambiar de sistema.'); return; }
if (UI.tour) UI.stopTour(true);
UI.select(null, { keepCam: true }); UI.openSection(null);
[Compare, Planner, TL].forEach(T => T.close ? T.close() : T.toggle && T.toggle(false)); Sci.closePop && Sci.closePop();
id === 'gargantua' ? this.enterGargantua() : this.enterSolar();
},
fade(on, ms) { const f = $('#sys-fade') || this.mkFade(); f.style.transitionDuration = ms + 'ms'; f.classList.toggle('on', on); },
mkFade() { const f = document.createElement('div'); f.id = 'sys-fade'; f.innerHTML = '<div class="sys-title"><b id="sys-t1"></b><span id="sys-t2"></span></div><button class="txt-btn sys-skip" id="sys-skip">Omitir</button>'; $('#app').appendChild(f); $('#sys-skip').addEventListener('click', () => this.skip && this.skip()); return f; },
apply(id) {
World.setSystem(id);
const app = $('#app'); Object.keys(SYSTEMS).forEach(k => app.classList.toggle('sys-' + k, k === id));
UI.SECTIONS = SECTIONS_BY_SYSTEM[id] || SECTIONS_SOLAR;
$('#rail').innerHTML = UI.SECTIONS.map(([sid, label]) => `<button class="rail-btn" data-sec="${sid}" aria-label="${label}" data-tip="${label}" data-tip-side="right">${ICON[sid] || ICON.herramientas}<span>${label}</span></button>`).join('');
const sb = $('#sys-btn'); if (sb) sb.querySelector('span').textContent = SYSTEMS[id].name;
},
enterGargantua() {
const G = SYSTEMS.gargantua, reduce = UI.reducedMotion, f = $('#sys-fade') || this.mkFade(), app = $('#app');
this.busy = true; SFX.play && SFX.play('select');
this.fade(true, reduce ? 250 : 900);
const T = [];
const done = () => { T.forEach(clearTimeout); this.busy = false; this.skip = null; Cam.introSY = 0; app.classList.remove('sys-intro'); f.classList.remove('on', 'titled'); World.bhReveal = 1;
Cam.fly = null; Cam.focus = World.sun; Cam.follow = true; Cam.target = [0, 0, 0]; Cam.dDist = Cam.dist = G.H * 34; Cam.el = 0.12; UI.hint('Sistema Gargantúa: elige un objeto en el menú o en el buscador'); };
T.push(setTimeout(() => {
this.apply('gargantua'); app.classList.add('sys-intro'); Cam.introSY = -0.33;     // Gargantúa en el tercio superior, título en el inferior
World.bhReveal = 0; Cam.fly = null; Cam.focus = World.sun; Cam.follow = true; Cam.target = [0, 0, 0]; Cam.az = 0.55; Cam.el = 0.05; Cam.dDist = Cam.dist = G.H * 150;
this.fade(false, reduce ? 300 : 1800);                                   // 2. aparecen las estrellas
if (reduce) { done(); return; }
const t0 = performance.now(), rev = () => { const k = clamp((performance.now() - t0 - 900) / 2600, 0, 1); World.bhReveal = k * k * (3 - 2 * k); if (k < 1 && this.busy) requestAnimationFrame(rev); };
requestAnimationFrame(rev);                                              // 3–4. Gargantúa y su disco se hacen visibles
T.push(setTimeout(() => Cam.travelTo([0, 0, 0], G.H * 34, 0.12), 1600));   // 5. aproximación lenta
T.push(setTimeout(() => { $('#sys-t1').textContent = 'SISTEMA GARGANTÚA'; $('#sys-t2').textContent = 'Interstellar'; f.classList.add('titled'); }, 2600));   // 6. título
T.push(setTimeout(() => f.classList.remove('titled'), 6200));
T.push(setTimeout(done, 7000));                                           // 7. controles
}, reduce ? 260 : 950));
this.skip = done;
},
enterSolar() {
const reduce = UI.reducedMotion; this.busy = true;
this.fade(true, reduce ? 250 : 800);
setTimeout(() => {
this.apply('solar'); Cam.fly = null; Cam.focus = World.byId.sol; Cam.follow = true; Cam.target = [0, 0, 0]; Cam.dDist = Cam.dist = ScaleState.get('overview'); Cam.el = 0.5;
this.fade(false, reduce ? 300 : 1100); this.busy = false; UI.toast('De vuelta en el Sistema Solar');
}, reduce ? 260 : 850);
},
};
const _secSys = UI.sectionHTML;
UI.sectionHTML = function (id) {
if (id === 'gargantua') {
const G = SYSTEMS.gargantua;
return `<p class="lead">${esc(G.sub)}. Explora el agujero negro supermasivo y los planetas que visita la misión Endurance.</p>
<div class="list">
<button class="cta-row" data-sys-free>${ICON.overview}<span><b>Explorar libremente</b><small>Navega por el sistema a tu ritmo</small></span></button>
<button class="cta-row ms-cta" data-ms-brief="endurance">${ICON.tour}<span><b>Modo misión Endurance</b><small>Gargantúa → Miller → Mann → maniobra gravitacional → Edmunds</small></span></button>
</div>
<h3>Objetos del sistema</h3>
<div class="list">${World.systems.gargantua.rb.map(rb => this.bodyRow(rb, rb.def.type === 'blackhole' ? 'Agujero negro supermasivo' : rb.def.info.fiction[0][1])).join('')}</div>
<p class="note sys-note">${esc(G.note)}</p><p class="note">${esc(G.orbitsNote)}</p>
<button class="cta-row" data-sys-go="solar">${ICON.overview}<span><b>Regresar al Sistema Solar</b><small>Volver a nuestro sistema planetario</small></span></button>`;
}
const h = _secSys.call(this, id);
if (id === 'explorar') return `<button class="cta-row sys-cta" data-sys-go="gargantua">${ICON.gargantua}<span><b>Explorar Sistema Gargantúa</b><small>Universo de Interstellar · sistema ficticio</small></span></button>` + h;
return h;
};
document.addEventListener('click', e => { const g = e.target.closest('[data-sys-go]'); if (g) Systems.go(g.dataset.sysGo); if (e.target.closest('#sys-btn')) Systems.openPicker(); });
const _renderInfo = UI.renderInfo;
UI.renderInfo = function (rb) {
if (rb.def.system !== 'gargantua') return _renderInfo.call(this, rb);
const d = rb.def, i = d.info, P = $('#info-body'), row = ([k, v]) => `<div class="row"><dt>${esc(k)}</dt><dd>${esc(v)}</dd></div>`;
P.innerHTML = `
<header class="info-head">
<span class="kind"><i class="dot" style="--c:${d.color}"></i>${esc(TYPE_LABEL[d.type])}${d.fiction || d.type === 'blackhole' ? ' · Interstellar' : ''}</span>
<h2 id="info-title">${esc(d.name)}</h2><p class="desc">${esc(i.desc)}</p>
<div class="info-actions"><button class="pill" id="info-travel">${ICON.travel}<span>Viajar</span></button><button class="pill" id="info-follow">${ICON.follow}<span>Seguir</span></button><button class="pill" id="info-close2">${ICON.zoomin}<span>Acercar</span></button></div>
</header>
<section class="gx-fic"><h3><span class="gx-tag fic">Universo de Interstellar</span></h3><dl>${i.fiction.map(row).join('')}</dl>${i.film ? `<p class="gx-film">${esc(i.film)}</p>` : ''}</section>
${i.dilation ? Dilation.html() : ''}
${i.science ? `<section class="gx-sci"><h3><span class="gx-tag sci">Ciencia real</span></h3>${i.science.map(([k, v]) => `<div class="gx-concept"><b>${esc(k)}</b><p>${esc(v)}</p></div>`).join('')}${i.render ? `<p class="note">${esc(i.render)}</p>` : ''}</section>` : ''}
<footer class="src"><p><b>Fuente:</b> ${esc(d.src)}</p><p class="sys-note">${esc(SYSTEMS.gargantua.note)}</p>${rb.def.orbit ? `<p class="note">${esc(SYSTEMS.gargantua.orbitsNote)}</p>` : ''}</footer>`;
P.scrollTop = 0;
$('#info-travel').addEventListener('click', () => { this.userActed(); Cam.travel(rb); });
$('#info-follow').addEventListener('click', () => { this.userActed(); if (Cam.focus === rb && Cam.follow) Cam.follow = false; else Cam.travel(rb, { keepAngles: true }); this.refreshCamButtons(); });
$('#info-close2').addEventListener('click', () => { this.userActed(); Cam.travel(rb, { close: true }); });
if (i.dilation) Dilation.start();
this.refreshCamButtons();
};
const _updateLive = UI.updateLive;
UI.updateLive = function (force) { if (this.sel && this.sel.def.system === 'gargantua') return; return _updateLive.call(this, force); };
const Dilation = {
RATIO: 7 * 365.25 * 24,                          // horas externas por hora local ≈ 61,362
t: 0, run: true, last: 0, timer: 0,
html() {
const ex = [[600, '10 minutos'], [3600, '1 hora'], [3 * 3600, '3 horas']];
return `<section class="gx-dil" id="info-dil"><h3>Dilatación temporal <span class="gx-tag fic">Universo de Interstellar</span></h3>
<p>En Miller, la intensa gravedad asociada a su proximidad a Gargantúa provoca una enorme diferencia entre el tiempo local y el tiempo experimentado lejos del planeta dentro del escenario de Interstellar.</p>
<p class="gx-ratio">1 hora en Miller ≈ 7 años externos <small>(factor ≈ ${fmt(MILLER_TIME.RATIO, 0)})</small></p>
<div class="gx-clocks"><div><span>Tiempo en Miller</span><b id="dil-local">00:00:00</b></div><div><span>Tiempo externo equivalente</span><b id="dil-ext">0 horas</b></div></div>
<p class="gx-ex">${ex.map(([s, t]) => `${t} → ${MillerFX.ext(s)}`).join(' · ')}</p>
<p class="note">El contador sigue el reloj de la ambientación sonora de Miller: cada tic llega cada ${fmt(MILLER_TIME.TICK_T, 4)} s locales (≈ ${fmt(MILLER_TIME.TICK_T * MILLER_TIME.RATIO / 3600, 1)} horas externas). La dilatación temporal gravitacional es real; esta relación concreta pertenece a la película.</p></section>`;
},
start() { /* sin reloj propio: lo actualiza MillerFX con el reloj de audio */ },
paint() {
const el = $('#dil-local'), ex = $('#dil-ext'); if (!el) return;
const s = Math.floor(this.t), p2 = n => String(n).padStart(2, '0');
el.textContent = `${p2(Math.floor(s / 3600))}:${p2(Math.floor(s / 60) % 60)}:${p2(s % 60)}`;
const days = this.t * this.RATIO / 86400; let y = Math.floor(days / 365.25), m = Math.round((days - y * 365.25) / 30.44); if (m === 12) { y++; m = 0; } const dd = Math.max(0, Math.round(days - y * 365.25 - m * 30.44));
ex.textContent = days < 1 ? fmt(days * 24, 1) + ' horas' : '≈ ' + [y ? y + (y === 1 ? ' año' : ' años') : '', m ? m + (m === 1 ? ' mes' : ' meses') : '', !y && dd ? dd + (dd === 1 ? ' día' : ' días') : ''].filter(Boolean).join(' y ');
},
};
(() => {
const b = document.createElement('button'); b.id = 'sys-btn'; b.className = 'sys-btn'; b.setAttribute('aria-haspopup', 'dialog'); b.dataset.tip = 'Cambiar de sistema';
b.innerHTML = ICON.systems + '<span>Sistema Solar</span>'; const fly = $('#fly-btn'); fly.parentNode.insertBefore(b, fly);
const w = document.createElement('div'); w.id = 'bh-warn'; w.className = 'ui'; w.setAttribute('role', 'status'); w.innerHTML = '<b>PROXIMIDAD A GARGANTÚA</b><span>Campo gravitacional extremo · la cámara no puede acercarse más al horizonte de sucesos</span>'; $('#app').appendChild(w);
const chip = document.createElement('div'); chip.id = 'sys-chip'; chip.className = 'ui'; chip.innerHTML = '<b>SISTEMA GARGANTÚA</b><span>Universo de Interstellar · sistema ficticio</span>'; $('#app').appendChild(chip);
})();
Systems.frame = function () {
if (World.system !== 'gargantua') return;
const bh = World.sun, d = V.len(V.sub(World.cam, bh.posS)) / bh.rS, near = d < 7 && !Systems.busy;
const w = $('#bh-warn'); if (w.classList.contains('on') !== near) w.classList.toggle('on', near);
World.bhBoost = clamp((9 - d) / 5, 0, 1) * 0.35;               // el disco se intensifica al acercarse
};
const MissionCam = {
pose() { return { target: Cam.target.slice(), dist: Cam.dist, az: Cam.az, el: Cam.el }; },
apply(p) { Cam.fly = null; Cam.focus = null; Cam.follow = false; Cam.auto = 0; Cam.target = p.target.slice(); Cam.dist = Cam.dDist = p.dist; Cam.az = p.az; Cam.el = p.el; },
blend(a, b, k) {
const e = k < 0.5 ? 16 * k ** 5 : 1 - Math.pow(-2 * k + 2, 5) / 2;          // arranque y llegada muy suaves
let daz = ((b.az - a.az) % (Math.PI * 2) + Math.PI * 3) % (Math.PI * 2) - Math.PI;
return { target: V.add(a.target, V.scale(V.sub(b.target, a.target), e)), dist: Math.exp(lerp(Math.log(a.dist), Math.log(b.dist), e)), az: a.az + daz * e, el: lerp(a.el, b.el, e) };
},
angles(dir) { const d = V.norm(dir); return { az: Math.atan2(d[0], d[2]), el: Math.asin(clamp(d[1], -1, 1)) }; },
};
const Missions = {
registry: {}, run: null, speed: 1,
register(def) { this.registry[def.id] = def; },
briefing(id) {
const M = this.registry[id]; if (!M) return;
UI.openSection(null); UI.select(null, { keepCam: true });
this.overlay(`<span class="ms-kicker">${esc(M.kicker)}</span><h2>${esc(M.title)}</h2><p class="ms-sub">${esc(M.sub)}</p><p class="ms-desc">${esc(M.desc)}</p>
<ol class="ms-route">${M.stages.filter(s => s.label).map(s => `<li>${esc(s.label)}</li>`).join('')}</ol><p class="note">${esc(M.note)}</p>
<div class="ms-btns"><button class="cta" data-ms="start" data-id="${id}">Iniciar misión</button><button class="cta-ghost" data-ms="close">Volver</button></div>`);
},
overlay(html) {
let o = $('#ms-ov'); if (!o) { o = document.createElement('section'); o.id = 'ms-ov'; o.setAttribute('role', 'dialog'); $('#app').appendChild(o); }
o.innerHTML = `<div class="ms-card">${html}</div>`; o.classList.add('open');
},
closeOverlay() { const o = $('#ms-ov'); if (o) o.classList.remove('open'); },
start(id) {
const M = this.registry[id]; if (!M) return;
this.closeOverlay();
if (this.run) this.exit({ restore: true, silent: true });
const saved = { cam: Object.assign(MissionCam.pose(), { focus: Cam.focus, follow: Cam.follow, fov: Cam.fov }), time: { jd: Time.jd, paused: Time.paused, idx: Time.idx, live: Time.live },
sel: UI.sel, mode: S.mode, section: UI.section };
if (S.mode === 'cine') UI.setMode('explore');
UI.select(null, { keepCam: true }); UI.openSection(null);
Time.paused = true;                                                       // las posiciones no avanzan durante la misión
this.run = { M, i: -1, t: 0, paused: false, saved, cine: false, hidden: false, shot: null, path: null, pathA: 0, ship: null };
$('#app').classList.add('in-mission'); this.hud(); this.go(0);
this.sfx('start');
},
exit(opts) {
opts = opts || {}; const R = this.run; if (!R) return;
const st = R.M.stages[R.i]; if (st && st.exit) st.exit(this.ctx);
this.run = null; this.ctx.cleanup(); if (MillerFX.active) MillerFX.exit();
$('#app').classList.remove('in-mission', 'mission-cine', 'mission-hide', 'ms-active');
const s = R.saved;
Time.jd = s.time.jd; Time.paused = s.time.paused; Time.idx = s.time.idx; Time.live = s.time.live; UI.refreshTime();
Cam.fov = s.cam.fov;
if (opts.keepCamera) { Cam.focus = opts.focus || null; Cam.follow = !!opts.focus; }
else { MissionCam.apply(s.cam); Cam.focus = s.cam.focus; Cam.follow = s.cam.follow; }
if (s.mode === 'cine' && !opts.keepCamera) UI.setMode('cine');
if (opts.select) UI.select(opts.select, { keepCam: true }); else if (s.sel && !opts.keepCamera && World.rb.includes(s.sel)) UI.select(s.sel, { keepCam: true });
['ms-hud', 'ms-ctl', 'ms-cap', 'ms-warn', 'ms-dil'].forEach(id => { const e = $('#' + id); if (e) e.remove(); });
this.closeOverlay();
if (!opts.silent) UI.toast(opts.msg || 'Has salido de la misión: vuelves al explorador libre');
},
go(i) {
const R = this.run; if (!R) return;
i = clamp(i, 0, R.M.stages.length - 1);
const prev = R.M.stages[R.i]; if (prev && prev.exit) prev.exit(this.ctx);
R.i = i; R.t = 0; R.shot = null; R.marks = {};
const st = R.M.stages[i]; this.caption(null);
if (st.enter) st.enter(this.ctx);
if (st.sfx) this.sfx(st.sfx);
this.paintHud();
},
next() { const R = this.run; if (R && R.i < R.M.stages.length - 1) this.go(R.i + 1); },
prev() { const R = this.run; if (R) this.go(Math.max(0, R.i - 1)); },
restart() { const R = this.run; if (R) { R.paused = false; this.go(0); } },
pause(v) {
const R = this.run; if (!R) return; R.paused = v == null ? !R.paused : v;
if (!R.paused && R.shot) { R.shot.from = MissionCam.pose(); R.shot.t0 = R.t; R.shot.dur = Math.min(R.shot.dur, 1.4); }   // reanudar con suavidad
this.paintHud();
},
toggleCine() { const R = this.run; if (!R) return; R.cine = !R.cine; $('#app').classList.toggle('mission-cine', R.cine); this.wake(); this.paintHud(); },
toggleHide() { const R = this.run; if (!R) return; R.hidden = !R.hidden; $('#app').classList.toggle('mission-hide', R.hidden); this.wake(); this.paintHud(); },
wake() { $('#app').classList.add('ms-active'); clearTimeout(this.wakeT); this.wakeT = setTimeout(() => { if (this.run && (this.run.cine || this.run.hidden) && !this.run.paused) $('#app').classList.remove('ms-active'); }, 3200); },
frame() {
const R = this.run; if (!R) return;
const now = performance.now(), dt = Math.min(0.1, (now - (this.lastT || now)) / 1000) * this.speed; this.lastT = now;
if (R.paused) return;
const st = R.M.stages[R.i];
R.t += dt;
if (R.shot) { const S = R.shot, k = clamp((R.t - S.t0) / S.dur, 0, 1), to = typeof S.to === 'function' ? S.to() : S.to; MissionCam.apply(MissionCam.blend(S.from, to, k)); }
if (st.update) st.update(this.ctx, R.t, dt);
if (st.dur && R.t >= st.dur && st.auto !== false) this.next();
if (R.pathA > 0 && !R.pathKeep) R.pathA = Math.max(0, R.pathA - dt * 0.4);
},
hud() {
const h = document.createElement('div'); h.id = 'ms-hud'; $('#app').appendChild(h);
const c = document.createElement('div'); c.id = 'ms-ctl'; c.setAttribute('role', 'toolbar'); c.setAttribute('aria-label', 'Controles de la misión'); $('#app').appendChild(c);
c.addEventListener('click', e => { const b = e.target.closest('[data-mc]'); if (!b) return; ({ prev: () => this.prev(), play: () => this.pause(), next: () => this.next(), restart: () => this.restart(), cine: () => this.toggleCine(), hide: () => this.toggleHide(), exit: () => this.exit({ restore: true }) })[b.dataset.mc](); this.wake(); });
const cap = document.createElement('div'); cap.id = 'ms-cap'; $('#app').appendChild(cap);
const w = document.createElement('div'); w.id = 'ms-warn'; $('#app').appendChild(w);
['pointermove', 'pointerdown', 'keydown'].forEach(t => addEventListener(t, () => { if (this.run) this.wake(); }, { passive: true }));
this.wake();
},
paintHud() {
const R = this.run; if (!R) return;
const st = R.M.stages, cur = st[R.i], steps = st.filter(s => s.label), idx = steps.indexOf(cur) >= 0 ? steps.indexOf(cur) : steps.length - 1;
$('#ms-hud').innerHTML = `<b>${esc(R.M.title)}</b><span>${idx + 1} / ${steps.length}</span><ol>${steps.map((s, k) => `<li class="${k < idx ? 'done' : k === idx ? 'on' : ''}">${esc(s.label)}</li>`).join('')}</ol>`;
const ic = (n, t, on) => `<button data-mc="${n}" class="${on ? 'on' : ''}" aria-label="${t}" data-tip="${t}">${MS_ICON[n]}</button>`;
$('#ms-ctl').innerHTML = ic('prev', 'Etapa anterior') + ic('play', R.paused ? 'Continuar' : 'Pausar').replace(MS_ICON.play, R.paused ? MS_ICON.resume : MS_ICON.play) + ic('next', 'Siguiente etapa') + ic('restart', 'Reiniciar misión')
+ '<i class="ms-sep"></i>' + ic('cine', 'Modo Cine', R.cine) + ic('hide', 'Ocultar interfaz', R.hidden) + ic('exit', 'Salir de la misión');
},
caption(title, sub, ms) {
const c = $('#ms-cap'); if (!c) return; clearTimeout(this.capT);
if (!title) { c.classList.remove('on'); return; }
c.innerHTML = `<b>${esc(title)}</b>${sub ? `<span>${esc(sub)}</span>` : ''}`; c.classList.add('on');
if (ms) this.capT = setTimeout(() => c.classList.remove('on'), ms / this.speed);
},
warn(t) { const w = $('#ms-warn'); if (!w) return; if (w.dataset.t === (t || '')) return; w.dataset.t = t || ''; w.textContent = t || ''; w.classList.toggle('on', !!t); },
sfx(n) {
if (!SFX.ready || !SFX.ready()) return; const T = (o) => SFX.tone('navigation', o), H = (o) => SFX.hiss('navigation', o);
try {
({ start: () => { T({ f: 110, f2: 82, dur: 2.4, gain: 0.05, lp: 700 }); T({ f: 165, f2: 123, dur: 2.4, gain: 0.03, lp: 900, delay: 0.15 }); },
whoosh: () => H({ f0: 900, dur: 1.6, gain: 0.03, q: 0.8 }),
deep: () => { T({ f: 52, f2: 41, dur: 5, gain: 0.07, lp: 300 }); H({ f0: 220, dur: 4, gain: 0.02, q: 0.6 }); },
assist: () => { T({ f: 55, f2: 110, dur: 7, gain: 0.06, lp: 500 }); H({ f0: 400, dur: 6, gain: 0.03, q: 0.7 }); },
arrive: () => { T({ f: 330, dur: 2.2, gain: 0.03, lp: 1600 }); T({ f: 494, dur: 2.2, gain: 0.025, lp: 1600, delay: 0.25 }); T({ f: 659, dur: 2.6, gain: 0.02, lp: 1800, delay: 0.5 }); },
cold: () => T({ f: 220, f2: 196, dur: 3, gain: 0.025, lp: 1200 }) })[n]?.();
} catch (e) { /* el audio es opcional */ }
},
};
const MS_ICON = {
prev: '<svg viewBox="0 0 24 24"><path d="M15 6l-6 6 6 6"/></svg>', next: '<svg viewBox="0 0 24 24"><path d="M9 6l6 6-6 6"/></svg>',
play: '<svg viewBox="0 0 24 24"><path d="M8.5 6v12M15.5 6v12"/></svg>', resume: '<svg viewBox="0 0 24 24"><path d="M8 5.5v13l10.5-6.5z"/></svg>',
restart: '<svg viewBox="0 0 24 24"><path d="M4.5 12a7.5 7.5 0 1 0 2.2-5.3M4.5 4.5v3.7h3.7"/></svg>',
cine: '<svg viewBox="0 0 24 24"><rect x="3" y="5" width="18" height="14" rx="2"/><path d="M7 5v14M17 5v14"/></svg>',
hide: '<svg viewBox="0 0 24 24"><path d="M3 12s3.5-6 9-6 9 6 9 6-3.5 6-9 6-9-6-9-6z"/><path d="M4 4l16 16"/></svg>',
exit: '<svg viewBox="0 0 24 24"><path d="M6 6l12 12M18 6L6 18"/></svg>',
};
Missions.ctx = {
get t() { return Missions.run ? Missions.run.t : 0; },
rb: id => World.byId[id],
shot(to, dur) { const R = Missions.run; R.shot = { from: MissionCam.pose(), to, t0: R.t, dur: Math.max(dur, 0.01) }; },
poseBehind(rb, bg, distK, daz, del) { const a = MissionCam.angles(V.sub(rb.posS, bg.posS)); return { target: rb.posS.slice(), dist: rb.rS * distK, az: a.az + (daz || 0), el: a.el + (del || 0) }; },
caption: (t, s, ms) => Missions.caption(t, s, ms), warn: t => Missions.warn(t), sfx: n => Missions.sfx(n),
once(key, fn) { const R = Missions.run; if (!R.marks[key]) { R.marks[key] = true; fn(); } },
setPath(pts) { const R = Missions.run; R.path = pts; R.pathA = 1; R.pathKeep = true; },
fadePath() { const R = Missions.run; if (R) R.pathKeep = false; },
ship(v) { const R = Missions.run; if (!v) { R.ship = null; return; } R.ship = R.ship || { pos: [0, 0, 0], vel: [0, 0, 1], size: 0.9, thr: 1 }; Object.assign(R.ship, v); },
grade(c) { $('#gl').style.filter = c || ''; },
dilation(on, factor) {
let d = $('#ms-dil');
if (!on) { if (d) d.remove(); return; }
if (!d) { d = document.createElement('div'); d.id = 'ms-dil'; $('#app').appendChild(d); }
d.dataset.f = factor; d.dataset.t = d.dataset.t || 0;
},
cleanup() {
World.bhBoost = 0; World.bhReveal = 1; Cam.introSY = 0; this.grade(''); Missions.warn('');
const d = $('#ms-dil'); if (d) d.remove();
},
};
Missions.draw = function (useP, t) {
const R = this.run; if (!R) return;
const cam = World.cam, P = World.P;
if (R.ship) {
const sh = R.ship; if (!sh.mesh) { PackedModels.request && PackedModels.request('ranger'); sh.mesh = CraftModels.mesh('ranger'); }
const m0 = CraftModels.mesh('ranger'); if (m0) sh.mesh = m0;
if (sh.mesh) {
const F = V.norm(sh.vel), U0 = Math.abs(F[1]) > 0.95 ? [1, 0, 0] : [0, 1, 0], X = V.norm(V.cross(U0, F)), U = V.cross(F, X);
const tex = !!sh.mesh.groups, pr = useP(tex ? P.craftTex : P.craft), rel = V.sub(sh.pos, cam);
GLX.state({ depthTest: true, depthWrite: true, cull: 'back', blend: 'none' });
const ep = V.add(rel, V.scale(F, -sh.size * 1.1));
GLX.setAll(pr, { u_model: M4.fromBasis(X, U, F, sh.size, rel), u_occ: new Float32Array(16), u_occN: 0, u_sunR: World.shadowSunR, u_ambient: S.ambient + 0.02, u_hover: 0,
u_eng: [...ep, sh.thr * 1.3], u_engCol: [0.55, 0.75, 1.0], u_engR: sh.size * 0.5 });
if (tex) World.drawGroups(sh.mesh, pr); else GLX.draw(sh.mesh, pr);
}
}
if (R.path && R.pathA > 0.01 && !R.cine) {
const pts = R.path, n = Math.min(pts.length, 200);
if (!this.pathMesh) { this._pP = new Float32Array(600); this._pF = new Float32Array(200); this.pathMesh = GLX.mesh({ a_pos: { data: this._pP, size: 3 }, a_frac: { data: this._pF, size: 1 } }, null, GLX.gl.LINE_STRIP, true); }
for (let i = 0; i < n; i++) { this._pP.set(pts[i], i * 3); this._pF[i] = i / (n - 1); }
GLX.update(this.pathMesh, 'a_pos', this._pP, n); GLX.update(this.pathMesh, 'a_frac', this._pF, n);
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const pr = useP(P.orbit);
GLX.setAll(pr, { u_model: M4.translate(V.sub([0, 0, 0], cam)), u_color: lin('#ffd9a0'), u_alpha: 0.85 * R.pathA, u_cur: (t * 0.1) % 1, u_fade: 0.65 });
GLX.draw(this.pathMesh, pr);
}
};
Missions.register((() => {
const H = () => SYSTEMS.gargantua.H, bh = () => World.byId.gargantua;
let path = null, peri = 0;
const buildAssist = () => {
const mann = World.byId.mann.posS, ed = World.byId.edmunds.posS;
const dIn = V.norm(V.scale(mann, -1)), dOut = V.norm(ed);
let ang = Math.acos(clamp(V.dot(dIn, dOut), -1, 1)); ang = clamp(ang, 0.7, 2.6);
const e = 1 / Math.sin(ang / 2), rp = H() * 12.5, p = rp * (1 + e);
const rh = V.norm(V.sub(dIn, dOut)), vh = V.norm(V.add(dIn, dOut));
const numax = Math.acos(-1 / e) - 0.04, pts = [];
for (let i = 0; i < 160; i++) {
const nu = -numax + 2 * numax * i / 159, r = p / (1 + e * Math.cos(nu));
if (r > H() * 70) continue;
const lift = H() * 1.8 * Math.cos(nu * 0.5);                                // ligeramente por encima del plano del disco
pts.push(V.add(V.add(V.scale(rh, r * Math.cos(nu)), V.scale(vh, r * Math.sin(nu))), [0, lift, 0]));
}
peri = pts.reduce((b, q, i) => V.len(q) < V.len(pts[b]) ? i : b, 0);
return pts;
};
const along = (pts, u) => { const f = clamp(u, 0, 1) * (pts.length - 1), i = Math.min(Math.floor(f), pts.length - 2), k = f - i; return { pos: V.add(pts[i], V.scale(V.sub(pts[i + 1], pts[i]), k)), vel: V.sub(pts[i + 1], pts[i]) }; };
return {
id: 'endurance', kicker: 'Sistema Gargantúa', title: 'MISIÓN ENDURANCE', sub: 'Gargantúa → Miller → Mann → maniobra gravitacional → Edmunds',
desc: 'Recreación cinematográfica e interactiva de la ruta de la Endurance por el sistema de Gargantúa. Puedes pausarla, avanzar o retroceder de etapa y salir en cualquier momento: el explorador libre se restaura tal como estaba.',
note: 'Recreación inspirada en Interstellar. Las trayectorias y los tiempos son cinematográficos; los conceptos de lente gravitacional, dilatación temporal y asistencia gravitacional son reales.',
stages: [
{ id: 'INTRO', label: 'Llegada', dur: 11, sfx: 'start',
enter(c) { World.bhReveal = 0; Cam.introSY = -0.33; MissionCam.apply({ target: [0, 0, 0], dist: H() * 170, az: 0.55, el: 0.05 }); c.shot({ target: [0, 0, 0], dist: H() * 52, az: 0.7, el: 0.1 }, 10); },
update(c, t) { const k = clamp((t - 0.8) / 4.5, 0, 1); World.bhReveal = k * k * (3 - 2 * k);
c.once('t1', () => c.caption('SISTEMA GARGANTÚA', 'Universo de Interstellar', 4200)); if (t > 5.6) c.once('t2', () => c.caption('MISIÓN ENDURANCE', 'Gargantúa → Miller → Mann → Edmunds', 4600)); },
exit() { World.bhReveal = 1; Cam.introSY = 0; } },
{ id: 'GARGANTUA', label: 'Gargantúa', dur: 15, sfx: 'deep',
enter(c) { c.shot({ target: [0, 0, 0], dist: H() * 15, az: 1.05, el: 0.09 }, 8); c.caption('GARGANTÚA', 'Agujero negro supermasivo', 6000); },
update(c, t) { World.bhBoost = clamp((t - 3) / 6, 0, 1) * 0.3; if (t > 9) c.once('orb', () => c.shot({ target: [0, 0, 0], dist: H() * 13, az: 1.05 + 0.55, el: 0.16 }, 6)); },
exit() { World.bhBoost = 0; } },
{ id: 'MILLER', label: 'Miller', dur: 24, sfx: 'whoosh',
enter(c) { const m = c.rb('miller'); c.shot(() => c.poseBehind(m, bh(), 7, 0.95, 0.16), 9); },
update(c, t) {
if (t > 8.5) c.once('cap', () => { c.caption('PLANETA DE MILLER', 'Dilatación temporal extrema · 1 hora local ≈ 7 años externos', 7000); MillerFX.enter({ mission: true }); });
if (t > 13) c.once('drift', () => c.shot(() => c.poseBehind(c.rb('miller'), bh(), 5.5, 1.2, 0.18), 10));
},
exit(c) { MillerFX.exit(); } },
{ id: 'MANN', label: 'Mann', dur: 26, sfx: 'whoosh',
enter(c) { const m = c.rb('miller'), mid = V.scale(V.add(m.posS, [0, 0, 0]), 0.55); c.shot({ target: mid, dist: H() * 26, az: Cam.az + 0.4, el: 0.22 }, 6); },
update(c, t) {
if (t > 6) c.once('go', () => { c.sfx('whoosh'); c.shot(() => c.poseBehind(c.rb('mann'), bh(), 7.5, -0.95, 0.14), 10); });
if (t > 14) c.once('cap', () => { c.caption('PLANETA DE MANN', 'Mundo helado · aislamiento', 7000); c.sfx('cold'); c.grade('saturate(0.72) hue-rotate(-8deg) brightness(0.96) contrast(1.05)'); });
if (t > 18) c.once('drift', () => c.shot(() => c.poseBehind(c.rb('mann'), bh(), 5.8, -1.25, 0.2), 8));
},
exit(c) { c.grade(''); } },
{ id: 'GRAVITY_ASSIST', label: 'Maniobra', dur: 30, sfx: 'assist',
enter(c) { path = buildAssist(); c.setPath(path); c.ship({ pos: path[0], vel: V.sub(path[1], path[0]), size: 0.9 });
c.shot(() => { const s = Missions.run.ship; return { target: s.pos, dist: 7, az: MissionCam.angles(V.sub(s.pos, [0, 0, 0])).az + 0.25, el: 0.12 }; }, 4);
c.caption('MANIOBRA GRAVITACIONAL', 'Recreación cinematográfica y educativa', 6000); },
update(c, t) {
const u = clamp((t - 1) / 26, 0, 1), q = along(path, u < 0.5 ? 0.5 * Math.pow(u * 2, 0.75) : 1 - 0.5 * Math.pow((1 - u) * 2, 0.75));   // más rápido cerca del perihelio
c.ship({ pos: q.pos, vel: q.vel });
const r = V.len(q.pos) / H();
World.bhBoost = clamp((26 - r) / 14, 0, 1) * 0.45; Cam.fov = 42 * Math.PI / 180 * (1 + clamp((22 - r) / 12, 0, 1) * 0.12);
c.warn(r < 14.5 ? 'MANIOBRA GRAVITACIONAL' : r < 24 ? 'CAMPO GRAVITACIONAL EXTREMO' : '');
if (t > 8) c.once('lat', () => { const n = V.norm(V.cross(V.sub(path[peri], path[0]), V.sub(path[path.length - 1], path[peri]))); const a = MissionCam.angles(V.add(n, [0, 0.25, 0])); c.shot({ target: path[peri], dist: H() * 30, az: a.az, el: clamp(a.el, -0.25, 0.25) }, 3); });
if (t > 15.5) c.once('post', () => c.shot(() => { const s = Missions.run.ship; return { target: s.pos, dist: 8, az: MissionCam.angles(s.vel).az, el: 0.1 }; }, 3));
if (t > 22.5) c.once('gen', () => { c.fadePath(); c.shot({ target: [0, 0, 0], dist: H() * 85, az: Cam.az + 0.3, el: 0.55 }, 5); });
},
exit(c) { c.warn(''); World.bhBoost = 0; Cam.fov = 42 * Math.PI / 180; Missions.run && (Missions.run.path = null); } },
{ id: 'EDMUNDS', label: 'Edmunds', dur: 26, sfx: null,
enter(c) { const e = c.rb('edmunds'); c.ship(null); c.shot(() => c.poseBehind(e, bh(), 32, 0.5, 0.18), 8); },
update(c, t) {
if (t > 7) c.once('appr', () => { c.caption('PLANETA DE EDMUNDS', 'Destino final de la misión', 7000); c.shot(() => c.poseBehind(c.rb('edmunds'), bh(), 6, 0.9, 0.2), 8); });
if (t > 15) c.once('orb', () => { c.sfx('arrive'); c.shot(() => c.poseBehind(c.rb('edmunds'), bh(), 4.6, 1.6, 0.12), 10); });
} },
{ id: 'COMPLETE', dur: 0, auto: false,
enter(c) { c.caption(null); Missions.overlay(`<span class="ms-kicker">Sistema Gargantúa</span><h2>MISIÓN ENDURANCE COMPLETADA</h2><p class="ms-sub">Destino final: Edmunds</p>
<div class="ms-btns ms-btns-v"><button class="cta" data-ms="explore">Explorar Edmunds</button><button class="cta-ghost" data-ms="again">Repetir misión</button><button class="cta-ghost" data-ms="system">Regresar al Sistema Gargantúa</button><button class="cta-ghost" data-ms="solar">Volver al Sistema Solar</button></div>`); } },
],
};
})());
document.addEventListener('click', e => {
const b = e.target.closest('[data-ms]'); if (!b) return; const a = b.dataset.ms;
if (a === 'start') Missions.start(b.dataset.id); else if (a === 'close') Missions.closeOverlay();
else if (a === 'explore') { const ed = World.byId.edmunds; Missions.exit({ keepCamera: true, focus: ed, select: ed, msg: 'Explorador libre: Planeta de Edmunds' }); }
else if (a === 'again') { Missions.closeOverlay(); Missions.restart(); }
else if (a === 'system') Missions.exit({ restore: true, msg: 'De vuelta al explorador del Sistema Gargantúa' });
else if (a === 'solar') { Missions.exit({ restore: true, silent: true }); Systems.go('solar'); }
});
document.addEventListener('click', e => { if (e.target.closest('[data-ms-brief]')) Missions.briefing(e.target.closest('[data-ms-brief]').dataset.msBrief); if (e.target.closest('[data-sys-free]')) UI.openSection(null); });
addEventListener('keydown', e => {
if (!Missions.run || e.target.closest && e.target.closest('input, textarea, select')) return;
if (e.key === 'Escape') { e.stopImmediatePropagation(); e.preventDefault(); Missions.exit({ restore: true }); }
else if (e.key === ' ') { e.stopImmediatePropagation(); e.preventDefault(); Missions.pause(); }
}, true);
const MILLER_TIME = {
RATIO: 7 * 365.25 * 24,              // horas externas por hora local (relación de Interstellar) ≈ 61,362
TICK_FIRST: 0.614, TICK_T: 1.2535,
};
const MillerFX = {
active: false, paused: false, src: null, nodes: null, t0: 0, pausedAt: 0, lastTick: -1, raf: 0, sched: 0, prevDuck: 1,
ctx() { try { return SFX.ensure ? SFX.ensure() : SFX.ctx; } catch (e) { return null; } },
level() { const A = Settings.state.audio; return (A.musicOn === false || A.muted) ? 0 : (A.music != null ? A.music : 0.6); },
clock() {
if (!this.active) return 0;
if (this.paused) return this.pausedAt;
if (this.src === 'track' && this.el) return this.el.currentTime;
const c = this.ctx(); return this.useAudio && c ? Math.max(0, c.currentTime - this.t0) : Math.max(0, performance.now() / 1000 - this.t0);
},
enter(opts) {
opts = opts || {};
if (this.active) { this.exit({ quick: true }); }                       // volver a Miller: la experiencia se reinicia sincronizada
this.active = true; this.paused = false; this.lastTick = -1; this.mission = !!opts.mission;
this.prevDuck = Music.duck; Music.setDuck(0);                           // fundido cruzado: la música principal baja
const c = this.ctx(), M = ASSET_MANIFEST.millerTrack;
if (M) this.startTrack(M); else this.startSynth(c);
this.hud(); this.loop();
},
exit(opts) {
opts = opts || {}; if (!this.active) return;
this.active = false; cancelAnimationFrame(this.raf); clearInterval(this.sched);
const c = this.ctx(), N = this.nodes, fade = opts.quick ? 0.25 : 1.6;
if (N && c) { const t = c.currentTime; N.out.gain.cancelScheduledValues(t); N.out.gain.setValueAtTime(N.out.gain.value, t); N.out.gain.linearRampToValueAtTime(0, t + fade);
setTimeout(() => { try { N.oscs.forEach(o => o.stop()); N.out.disconnect(); } catch (e) {} }, fade * 1000 + 100); }
if (this.el) { const el = this.el; const g = this.trackGain; if (g && c) g.gain.linearRampToValueAtTime(0, c.currentTime + fade); setTimeout(() => { el.pause(); }, fade * 1000); }
this.nodes = null;
Music.setDuck(this.prevDuck || 1);                                      // la música principal vuelve
const h = $('#mx-hud'); if (h) { h.classList.remove('on'); setTimeout(() => { if (!this.active) h.remove(); }, 700); }
},
pause(v) {
if (!this.active || v === this.paused) return;
const c = this.ctx();
if (v) { this.pausedAt = this.clock(); this.paused = true; clearInterval(this.sched); if (this.nodes && c) { this.nodes.out.gain.setTargetAtTime(0, c.currentTime, 0.08); this.tickBus.gain.setValueAtTime(0, c.currentTime); } if (this.el) this.el.pause(); this.tickUI(); }
else { this.paused = false; this.t0 = (this.useAudio && c ? c.currentTime : performance.now() / 1000) - this.pausedAt;
if (this.nodes && c) { this.nodes.out.gain.setTargetAtTime(this.level() * 0.55, c.currentTime, 0.15);
const nb = c.createGain(); nb.gain.value = 1; nb.connect(this.nodes.out); this.tickBus = nb;      // los tics ya programados quedan en el bus silenciado
this.nextTick = Math.ceil((this.pausedAt - MILLER_TIME.TICK_FIRST) / MILLER_TIME.TICK_T); }
if (this.el) this.el.play().catch(() => {}); else this.schedule(); }
},
startSynth(c) {
this.src = 'synth';
this.useAudio = !!(c && c.state === 'running');
if (!this.useAudio) { this.t0 = performance.now() / 1000; return; }       // audio aún no permitido: reloj del navegador, igual de continuo
const out = c.createGain(); out.gain.value = 0; out.connect(SFX.bus && SFX.bus.music ? SFX.bus.music : c.destination);
const lp = c.createBiquadFilter(); lp.type = 'lowpass'; lp.frequency.value = 900; lp.connect(out);
const oscs = [];
[[110, 'sine', 0.22], [164.81, 'sine', 0.13], [220.5, 'triangle', 0.06], [329.6, 'sine', 0.035]].forEach(([f, type, g]) => {
const o = c.createOscillator(), gn = c.createGain(), lfo = c.createOscillator(), lg = c.createGain();
o.type = type; o.frequency.value = f; gn.gain.value = g; lfo.frequency.value = 0.05 + Math.random() * 0.05; lg.gain.value = g * 0.45;
lfo.connect(lg); lg.connect(gn.gain); o.connect(gn); gn.connect(lp); o.start(); lfo.start(); oscs.push(o, lfo);
});
this.tickBus = c.createGain(); this.tickBus.gain.value = 1; this.tickBus.connect(out);
this.nodes = { out, oscs };
this.t0 = c.currentTime + 0.05;
out.gain.setValueAtTime(0, c.currentTime); out.gain.linearRampToValueAtTime(this.level() * 0.55, c.currentTime + 2.5);   // entrada suave
this.nextTick = 0; this.schedule();
},
schedule() {
clearInterval(this.sched);
const c = this.ctx(); if (!c || !this.nodes) return;
const run = () => {
if (!this.active || this.paused) return;
this.tickUI();                                                      // la interfaz no depende de los FPS del dibujo 3D
const local = this.clock(), horizon = local + 6;                    // antelación amplia: inmune a bloqueos largos del navegador
let k = Math.max(this.nextTick, Math.ceil((local - MILLER_TIME.TICK_FIRST) / MILLER_TIME.TICK_T));
for (; MILLER_TIME.TICK_FIRST + k * MILLER_TIME.TICK_T < horizon; k++) {
const at = this.t0 + MILLER_TIME.TICK_FIRST + k * MILLER_TIME.TICK_T; if (at < c.currentTime) continue;
const o = c.createOscillator(), g = c.createGain(), bp = c.createBiquadFilter();
o.type = 'square'; o.frequency.value = 2100; bp.type = 'bandpass'; bp.frequency.value = 3200; bp.Q.value = 6;
g.gain.setValueAtTime(0, at); g.gain.linearRampToValueAtTime(0.09, at + 0.0015); g.gain.exponentialRampToValueAtTime(0.0008, at + 0.035);
o.connect(bp); bp.connect(g); g.connect(this.tickBus); o.start(at); o.stop(at + 0.05);
(this.sched_log = this.sched_log || []).push(+(at - this.t0).toFixed(4)); if (this.sched_log.length > 32) this.sched_log.shift();   // diagnóstico
}
this.nextTick = k;
};
run(); this.sched = setInterval(run, 90);
},
startTrack(url) {
this.src = 'track';
const el = this.el || (this.el = new Audio()); el.src = url; el.preload = 'auto'; el.currentTime = 0;
const c = this.ctx();
if (c && !this.trackNode) { try { this.trackNode = c.createMediaElementSource(el); this.trackGain = c.createGain(); this.trackNode.connect(this.trackGain); this.trackGain.connect(SFX.bus.music); } catch (e) {} }
if (this.trackGain && c) { this.trackGain.gain.setValueAtTime(0, c.currentTime); this.trackGain.gain.linearRampToValueAtTime(this.level(), c.currentTime + 2.5); }
if (el.readyState < 3) this.loading(true); el.oncanplay = () => this.loading(false);
el.play().catch(() => {});
},
loading(v) { const h = $('#mx-hud'); if (h) h.classList.toggle('loading', v); },
hud() {
let h = $('#mx-hud');
if (!h) { h = document.createElement('div'); h.id = 'mx-hud'; h.setAttribute('role', 'status'); $('#app').appendChild(h); }
h.innerHTML = `<span class="mx-k">Dilatación temporal</span>
<div class="mx-row"><span>Miller</span><b id="mx-loc">00:00:00</b></div>
<div class="mx-row mx-ext"><span>Exterior</span><b id="mx-ext">0 horas</b></div>
<div class="mx-line"><i id="mx-dot"></i></div>
<small class="mx-rel">1 hora ≈ 7 años · relación de Interstellar</small>
<span class="mx-load">Preparando audio de Miller…</span>
<div class="mx-min"><b>MILLER</b><span>1 hora ≈ 7 años</span></div>`;
requestAnimationFrame(() => h.classList.add('on'));
},
ext(localS) {
const days = localS * MILLER_TIME.RATIO / 86400;
if (days < 1) return '≈ ' + fmt(days * 24, 1) + ' horas';
if (days < 60) return '≈ ' + fmt(days, 1) + ' días';
let y = Math.floor(days / 365.25), m = Math.round((days - y * 365.25) / 30.44); if (m === 12) { y++; m = 0; }
return '≈ ' + (y ? y + (y === 1 ? ' año' : ' años') + (m ? ' y ' : '') : '') + (m || !y ? m + (m === 1 ? ' mes' : ' meses') : '');
},
paint(local) {
const s = Math.floor(local), p2 = n => String(n).padStart(2, '0'), L = `${p2(Math.floor(s / 3600))}:${p2(Math.floor(s / 60) % 60)}:${p2(s % 60)}`, E = this.ext(local);
['mx-loc', 'dil-local'].forEach(id => { const e = $('#' + id); if (e && e.textContent !== L) e.textContent = L; });
['mx-ext', 'dil-ext'].forEach(id => { const e = $('#' + id); if (e && e.textContent !== E) e.textContent = E; });
const ph = ((local - MILLER_TIME.TICK_FIRST) % MILLER_TIME.TICK_T + MILLER_TIME.TICK_T) % MILLER_TIME.TICK_T / MILLER_TIME.TICK_T;
const d = $('#mx-dot'); if (d) d.style.transform = `translateX(${(ph * 100).toFixed(1)}%)`;
},
tickUI() {
if (!this.active) return;
const local = this.clock(), k = Math.floor((local - MILLER_TIME.TICK_FIRST) / MILLER_TIME.TICK_T);
this.paint(local);
if (k > this.lastTick && local >= MILLER_TIME.TICK_FIRST) {         // tic: pulsación sutil del tiempo externo
this.lastTick = k;
['mx-hud', 'info-dil'].forEach(id => { const e = $('#' + id); if (e) { e.classList.remove('tick'); void e.offsetWidth; e.classList.add('tick'); } });
}
},
loop() {
const step = () => { if (!this.active) return; this.tickUI(); this.raf = requestAnimationFrame(step); };
this.raf = requestAnimationFrame(step);
if (!this.nodes) { clearInterval(this.sched); this.sched = setInterval(() => this.tickUI(), 90); }   // sin audio: misma actualización periódica
},
};
const _selectMX = UI.select;
UI.select = function (rb, opts) {
const r = _selectMX.call(this, rb, opts);
if (Missions.run) return r;                                          // en la misión, la etapa MILLER controla la experiencia
if (rb && rb.id === 'miller') { if (!MillerFX.active) MillerFX.enter(); }
else if (MillerFX.active) MillerFX.exit();
return r;
};
const _msPause = Missions.pause.bind(Missions);
Missions.pause = function (v) { _msPause(v); if (MillerFX.active && Missions.run) MillerFX.pause(Missions.run.paused); };
const App = {
last: 0, scaleRatio: null, ready: false,
qualityParams() {
const g = Settings.state.graphics, d = Math.min(devicePixelRatio || 1, 2);
const lv = { low: 3, medium: 4 }[g.preset] || 5;
return { dpr: clamp(d * g.scale, 0.35, 2.5), samples: g.aa ? 4 : 0, bloomLevels: lv };
},
applyQuality(rebuild) {
const p = this.qualityParams();
GLX.samples = p.samples; GLX.bloomLevels = p.bloomLevels; this.dpr = p.dpr;
if (rebuild) { GLX.W = 0; this.resize(); }
},
resize() { GLX.resize(innerWidth, innerHeight, this.dpr); },
async start() {
const cv = $('#gl');
try { GLX.init(cv); }
catch (e) {
$('#load-msg').textContent = 'Este simulador necesita WebGL2. Actualiza el navegador o activa la aceleración por hardware en su configuración.';
$('#app').classList.add('failed'); return;
}
document.querySelectorAll('img[data-logo]').forEach(i => { i.src = document.querySelector(i.dataset.logo === 'big' ? '.intro-logo img' : '.about-logo').src; });
SettingsApply.boot();
I18N.set(Settings.state.language, true);
this.applyQuality(false); this.resize();
addEventListener('resize', () => this.resize());
Time.startState();
const steps = 9; let n = 0;
try {
await World.build(msg => { n++; $('#load-msg').textContent = msg; $('#load-bar').style.transform = `scaleX(${n / steps})`; });
} catch (e) {
console.error(e);
$('#load-msg').textContent = 'No se pudo preparar la escena 3D (' + e.message + '). Prueba con otro navegador o reduce la calidad.';
$('#app').classList.add('failed'); return;
}
UI.init(); Input.init(cv);
Cam.dist = Cam.dDist = ScaleState.get('overview') * 0.78; Cam.el = 0.34; Cam.az = 0.4; Cam.auto = UI.reducedMotion ? 0 : 0.03;
$('#load-msg').textContent = 'Listo';
$('#app').classList.add('ready');
const st = $('#start'); st.disabled = false; st.focus({ preventScroll: true });
st.addEventListener('click', () => this.begin());
const sg = $('#start-tours'); sg.disabled = false; sg.addEventListener('click', () => UI.toggleTours(true));
const sf = $('#start-fly'); sf.disabled = false; sf.addEventListener('click', () => { this.begin(); Flight.openHangar(); });
this.ready = true;
if (Settings.state.general.startup === 'explore' && location.hash !== '#creditos') this.begin();
window.Solaris = { World, Cam, UI, Time, S, select: id => UI.select(World.byId[id], { fly: true }), App, Flight, Music, SFX, Settings, Keys, Gfx, GLX, I18N, BODY, PackedModels, Assets, Compare, TL, Sci, ScaleState, Planner, Flight, Systems, Dilation, Missions, Cam: Cam, MillerFX, Music };
requestAnimationFrame(t => this.loop(t));
},
home() {
if (Flight.on) Flight.exit();
UI.stopTour();
if (S.mode !== 'explore') UI.setMode('explore');
UI.select(null); UI.openSection(null); UI.toggleHelp(false);
UI.started = false;
$('#app').classList.remove('started'); $('#app').classList.add('booting');
Time.startState(); UI.refreshTime();
Cam.travelTo([0, 0, 0], ScaleState.get('overview') * 0.78, 0.34);
Cam.auto = UI.reducedMotion ? 0 : 0.03;
setTimeout(() => $('#start').focus({ preventScroll: true }), 50);
},
begin(opts) {
opts = opts || {};
$('#app').classList.remove('booting'); $('#app').classList.add('started');
UI.started = true; Cam.auto = 0;
Time.startState(); UI.refreshTime();
if (!opts.quiet) {
Cam.travelTo([0, 0, 0], ScaleState.get('overview'), 0.42);
setTimeout(() => UI.hint(I18N.t(matchMedia('(pointer: coarse)').matches ? 'Toca cualquier objeto o búscalo por su nombre' : 'Haz clic en cualquier objeto o búscalo por su nombre'), 4200), 900);
}
$('#gl').focus({ preventScroll: true });
},
prepare() {
const W = innerWidth, H = innerHeight;
for (const rb of World.rb) {
const p = World.project(rb.posS, rb.proj);
p.rpx = p.on ? rb.rS * World.pxPerUnit(p.w) : 0;
const near = V.dist(rb.posS, Cam.pos) < rb.rS * 3;
const onScreen = p.on && p.x > -p.rpx - 40 && p.x < W + p.rpx + 40 && p.y > -p.rpx - 40 && p.y < H + p.rpx + 40;
const secondary = /^(moon|asteroid|comet|tno|craft)$/.test(rb.def.type);
rb.drawOn = (onScreen && p.rpx > (secondary ? Gfx.distPx() : 0.3)) || near;
if (rb.def.type === 'moon' && !S.layers.moons) rb.drawOn = false;
if (rb.isCraft && (rb.hidden || !S.layers.craft)) { rb.drawOn = false; p.on = false; p.rpx = 0; }
}
const sp = World.sun.proj;
World.starDim = 1 - 0.55 * (sp.on ? smoothstep(8, 220, sp.rpx) : 0);
},
advance(sec) {
this.hold = true;
const dt = 0.05;
if (Flight.on) { for (let t = 0; t < sec; t += dt) { Flight.simulate(dt); Flight.camera(dt); } Flight.forceHud = true; Flight.renderFrame(performance.now() / 1000); return; }
for (let t = 0; t < sec; t += dt) { Time.step(dt); World.update(Time.jd, dt, Time.speed); this.scaleCam(); Cam.update(dt); UI.tickTour(dt); this.prepare(); }
World.texTick(); PackedModels.tick();
World.render(performance.now() / 1000); World.frame++;
if (UI.started) { UI.updateLabels(); UI.updateLive(true); UI.tickTime(); UI.drawMinimap(); }
},
adapt(dt) {
this.ema = this.ema == null ? dt : this.ema * 0.93 + dt * 0.07;
this.adaptT = (this.adaptT || 0) + dt;
SettingsApply.autoCheck(dt);
if (!Settings.state.graphics.auto) return;     // la resolución adaptativa solo actúa con la calidad automática
if (this.adaptT < 1.5 || document.hidden) return;
this.adaptT = 0;
const max = this.qualityParams().dpr, min = Math.min(max, 0.75);
let d = this.dpr;
if (this.ema > 1 / 42 && d > min) d = Math.max(min, d * 0.85);
else if (this.ema < 1 / 57 && d < max) d = Math.min(max, d * 1.12);
if (Math.abs(d - this.dpr) > 0.01) { this.dpr = d; this.resize(); }
},
scaleCam() {
if (this.scaleRatio != null && !Cam.fly) {
const ref = Cam.focus ? Cam.focus.rS : ScaleState.get('overview');
Cam.dist = Cam.dDist = this.scaleRatio * ref;
if (ScaleState.t >= 1) this.scaleRatio = null;
}
},
loop(t) {
requestAnimationFrame(tt => this.loop(tt));
if (this.hold) return;
const cap = Settings.state.graphics.fps;
if (cap && this.last && t - this.last < 1000 / cap - 1.5) return;
const dt = this.last ? Math.min(0.1, (t - this.last) / 1000) : 0.016;
this.tickT = (this.tickT || 0) + dt;
if (this.tickT > 0.25 && World.texState) { this.tickT = 0; World.texTick(); PackedModels.tick(); }
this.last = t;
this.adapt(dt);
if (Flight.on) { Flight.frame(dt, t / 1000); return; }
Time.step(dt);
World.update(Time.jd, dt, Time.speed);
this.scaleCam();
Cam.update(dt);
this.prepare();
World.render(t / 1000);
World.frame++;
if (UI.started) UI.frame(dt);
},
};
window.addEventListener('DOMContentLoaded', () => App.start());
const SHIPS = [
{ id: 'colibri', name: 'Colibrí', type: 'Nave exploradora ligera', len: 22, maxSpeed: 300, accel: 60, angAcc: 110, maxAng: 80, mass: 18,
prop: 'Motor magnetoplasmático de alto empuje (concepto)', eng: [0.45, 0.68, 1.0],
desc: 'Pequeña, ágil y fácil de controlar. Ideal para aprender a pilotar y recorrer el entorno de un planeta o una luna.' },
{ id: 'ranger', name: 'Ranger One', type: 'Nave de exploración planetaria', len: 24, maxSpeed: 120, accel: 35, angAcc: 60, maxAng: 40, mass: 32,
prop: 'Motores de fusión compactos con propulsores vectoriales (concepto)', eng: [0.55, 0.75, 1.0],
info: { franquicia: 'Interstellar', clase: 'Lanzadera de exploración', funcion: 'Enlace entre la nave nodriza, las órbitas y la superficie de los planetas',
tripulacion: 'Reducida', afiliacion: 'Misión Endurance', destacadas: ['Capaz de operar en órbita y en atmósferas', 'Gran precisión en maniobras de acoplamiento'] },
brand: { src: 'assets/images/interstellar.png', alt: 'Interstellar', w: 640, h: 154, credit: { name: 'Logotipo de Interstellar', author: 'Titulares de los derechos de la película Interstellar' } },
desc: 'Una lanzadera de exploración inspirada en las naves de Interstellar, diseñada para operar entre planetas, estaciones orbitales y lunas con gran precisión. Combina agilidad, autonomía y capacidad de maniobra para afrontar desde operaciones orbitales hasta misiones de exploración en el espacio profundo.' },
{ id: 'rocinante', name: 'Rocinante', type: 'Corbeta de combate y patrulla', len: 46, maxSpeed: 3000, accel: 49, angAcc: 55, maxAng: 38, mass: 600, massEst: true, needsPack: true,
prop: 'Impulsor de fusión Epstein y propulsores de maniobra', eng: [0.72, 0.82, 1.0],
desc: 'Nave de alta maniobrabilidad adaptada para misiones de patrulla, combate, escolta y exploración táctica. Su diseño robusto y su capacidad operativa la convierten en una de las naves más reconocibles del universo de The Expanse.',
info: { franquicia: 'The Expanse', clase: 'Corbeta (fragata ligera) de origen marciano', funcion: 'Patrulla, combate, escolta y exploración táctica',
tripulacion: 'Reducida: en la serie la opera un equipo de entre cuatro y seis personas', velocidad: 'Sin viaje más rápido que la luz; aceleraciones sostenidas de varios g',
capacidad: 'Largo alcance dentro del Sistema Solar y maniobras de combate exigentes', armamento: 'Cañones de defensa puntual (PDC), lanzadores de torpedos y, en etapas posteriores, un cañón de riel en la quilla',
afiliacion: 'Construida para la armada de la República Congresional Marciana; después, nave independiente de su tripulación',
destacadas: ['Ágil y potente pese a su casco robusto', 'Pensada para aceleraciones altas y maniobras de combate', 'Una de las naves más reconocibles de la saga'] },
brand: { src: 'assets/images/logos/the-expanse.png', alt: 'The Expanse', w: 793, h: 232, credit: { name: 'Logotipo de The Expanse', author: 'Titulares de los derechos de la serie The Expanse' } } },
{ id: 'orville', name: 'ECV-197 Orville', type: 'Nave de exploración', len: 300, lenEst: true, maxSpeed: 30000, accel: 60, angAcc: 14, maxAng: 9, mass: 120000, massEst: true, needsPack: true,
prop: 'Impulso cuántico y motores subluz (en SOLARIS solo se simula el vuelo subluz)', eng: [0.45, 0.66, 1.0],
desc: 'Nave de exploración y servicio interestelar diseñada para misiones prolongadas, observación, contacto y apoyo operativo. Su perfil combina funciones científicas, diplomáticas y tácticas dentro del universo de The Orville.',
info: { franquicia: 'The Orville', clase: 'Nave de exploración de nivel medio de la Unión Planetaria', funcion: 'Exploración, observación, contacto y apoyo operativo',
tripulacion: 'Numerosa: cientos de tripulantes de distintas especies', velocidad: 'Viajes interestelares con impulso cuántico; en SOLARIS vuela solo a velocidades inferiores a la luz',
capacidad: 'Misiones prolongadas lejos de cualquier base', armamento: 'Defensivo, con armas de energía, torpedos y escudos', afiliacion: 'Unión Planetaria',
destacadas: ['Diseñada para viajes prolongados con tripulación numerosa', 'Multipropósito: ciencia, diplomacia y apoyo táctico', 'Perfil institucional, más de exploración que militar'] },
brand: { src: 'assets/images/logos/the-orville.svg', alt: 'The Orville', w: 1000, h: 140, credit: { name: 'Logotipo de The Orville', author: 'Titulares de los derechos de la serie The Orville' } } },
{ id: 'atlas', name: 'Atlas', type: 'Transbordador', len: 37, maxSpeed: 30, accel: 29, angAcc: 45, maxAng: 30, mass: 105,
prop: 'Motores químicos de metano y oxígeno líquidos', eng: [1.0, 0.6, 0.22],
desc: 'Prestaciones cercanas a la tecnología actual: perfecto para órbitas terrestres y viajes a la Luna, pero lento para ir más lejos.' },
{ id: 'odisea', name: 'Odisea', type: 'Nave de largo alcance', len: 140, maxSpeed: 8000, accel: 40, angAcc: 8, maxAng: 6, mass: 2400,
prop: 'Propulsión por fusión nuclear (concepto teórico)', eng: [0.78, 0.72, 1.0],
desc: 'Una nave de fusión capaz de cruzar el Sistema Solar en días. Pesada y lenta de maniobrar, pero incansable.' },
];
for (let i = SHIPS.length - 1; i >= 0; i--) if (SHIPS[i].needsPack && !MODEL_PACKS[SHIPS[i].id]) SHIPS.splice(i, 1);
const WARPS = [1, 10, 100, 1000, 10000, 100000];
const DETENTS = [-0.25, 0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00];
const DETENT_NAMES = { '-0.25': 'Inverso', '0': 'Motores en reposo', '0.1': 'Maniobra mínima', '0.25': 'Empuje bajo', '0.5': 'Crucero', '0.75': 'Empuje alto', '0.9': 'Potencia máxima continua', '1': 'Empuje máximo' };
const SPOOL = { colibri: [0.55, 1.1], ranger: [0.75, 0.95], atlas: [0.8, 0.9], odisea: [2.4, 0.38] };
function smoothDamp(cur, target, ref, smoothTime, maxSpeed, dt) {
if (dt <= 0) return cur;
smoothTime = Math.max(1e-4, smoothTime);
const omega = 2 / smoothTime, x = omega * dt, ex = 1 / (1 + x + 0.48 * x * x + 0.235 * x * x * x);
let change = cur - target; const orig = target, maxC = maxSpeed * smoothTime;
change = clamp(change, -maxC, maxC); target = cur - change;
const temp = (ref.v + omega * change) * dt;
ref.v = (ref.v - omega * temp) * ex;
let out = target + (change + temp) * ex;
if ((orig - cur > 0) === (out > orig)) { out = orig; ref.v = (out - orig) / dt; }
return out;
}
const CAMFX = {
cockpit: { fov: 64, fovUp: 4, back: 0.006, shake: 0.0035, steady: 0.0004 },
chase: { fov: 60, fovUp: 8, back: 0.12, shake: 0.006, steady: 0.0006 },
free: { fov: 55, fovUp: 3, back: 0.04, shake: 0.004, steady: 0.0004 },
cine: { fov: 47, fovUp: 10, back: 0.16, shake: 0.008, steady: 0.0008 },
front: { fov: 55, fovUp: 4, back: 0.06, shake: 0.004, steady: 0.0004 },
side: { fov: 55, fovUp: 4, back: 0.05, shake: 0.004, steady: 0.0004 },
};
const C_KMS = 299792.458;
const CAMS = [['cockpit', 'Cabina'], ['chase', 'Exterior'], ['free', 'Libre'], ['cine', 'Cinematográfica'], ['front', 'Frontal'], ['side', 'Lateral']];
const STARTS = [['tierra', 'Órbita terrestre'], ['luna', 'Órbita lunar'], ['marte', 'Órbita de Marte'], ['jupiter', 'Órbita de Júpiter'], ['saturno', 'Órbita de Saturno'], ['pluton', 'Órbita de Plutón']];
(() => {
const mirror = pts => pts.map(q => [-q[0], q[1]]);
const W = [0.8, 0.82, 0.86, 0.5], CAN = [0.04, 0.07, 0.11, 1.4], ACC = [0.92, 0.5, 0.12, 0.6], DK = [0.09, 0.09, 0.11, 0.9], NZ = [0.28, 0.28, 0.3, 0.7];
CraftModels.define('colibri', (B, C) => {
B.cyl([0, 0, 0], [0, 0, 1], 1.5, 13, 10, W, C.dark);
B.dish([0, 0, 6.5], [0, 0, -1], 1.5, 4.5, 10, W, W);
B.box([0, 1.15, 3.6], [1.5, 0.9, 3.4], CAN);
const wing = [[1.4, 2.5], [7.5, -3.5], [7.5, -5.2], [1.4, -4.5]];
B.poly([0, -0.3, -1.5], [1, 0, 0], [0, 0, 1], wing, W, C.dark, 0.25);
B.poly([0, -0.3, -1.5], [1, 0, 0], [0, 0, 1], mirror(wing), W, C.dark, 0.25);
[-1, 1].forEach(s => B.box([s * 7.5, -0.3, -5.6], [0.3, 0.9, 2.4], ACC));
B.poly([0, 1.0, -5.2], [0, 1, 0], [0, 0, 1], [[0, 1.6], [2.4, -0.4], [2.4, -1.4], [0, -1.6]], ACC, ACC, 0.18);
B.cyl([0, 0, -7.2], [0, 0, 1], 1.15, 1.6, 10, C.dark);
B.dish([0, 0, -9.4], [0, 0, -1], 1.05, 1.4, 12, C.dark, NZ);
});
CraftModels.define('ranger', (B, C) => {
B.box([0, 0, 0], [7, 2.6, 18], W); B.dish([0, 0, 9], [0, 0, -1], 2.2, 3, 10, W, W);
B.box([0, 1.4, 4], [3, 0.8, 4], CAN);
[-1, 1].forEach(s => B.poly([0, -0.5, -3], [1, 0, 0], [0, 0, 1], (s > 0 ? [[3.5, 3], [6.7, -4], [6.7, -7], [3.5, -7]] : [[-3.5, 3], [-6.7, -4], [-6.7, -7], [-3.5, -7]]), W, C.dark, 0.3));
[[-1, 0.3], [0, -1.2], [1, 0.3]].forEach(([x, y]) => B.dish([x, y, -12], [0, 0, -1], 0.7, 1.2, 10, C.dark, NZ));
});
CraftModels.define('atlas', (B, C) => {
B.box([0, 0, 0], [4.6, 4.2, 30], W);
B.box([0, -2.15, 0], [4.7, 0.3, 30], DK);
B.dish([0, 0, 15], [0, 0, -1], 2.3, 5, 12, W, DK);
B.box([0, 1.95, 12.6], [3.0, 0.6, 2.2], CAN);
const wing = [[2.3, 8], [12, -10], [12, -13], [2.3, -13]];
B.poly([0, -1.6, -2], [1, 0, 0], [0, 0, 1], wing, W, DK, 0.5);
B.poly([0, -1.6, -2], [1, 0, 0], [0, 0, 1], mirror(wing), W, DK, 0.5);
B.poly([0, 2.1, -11], [0, 1, 0], [0, 0, 1], [[0, 0], [7.5, -3.5], [7.5, -6], [0, -4]], W, W, 0.4);
[[0, 1.2], [1.4, -0.8], [-1.4, -0.8]].forEach(([x, y]) => B.dish([x, y, -16.8], [0, 0, -1], 0.95, 1.8, 12, C.dark, NZ));
});
CraftModels.define('odisea', (B, C) => {
B.box([0, 0, -5], [2, 2, 112], C.silver);
B.cyl([0, 0, 52], [0, 0, 1], 4.5, 14, 16, W, C.dark);
B.dish([0, 0, 59], [0, 0, -1], 4.5, 6, 16, W, W);
for (let k = 0; k < 20; k++) {
const a = TAU * k / 20, r = [Math.cos(a), Math.sin(a), 0], t = [-Math.sin(a), Math.cos(a), 0];
B.box([r[0] * 16, r[1] * 16, 34], [5.4, 3.2, 4], W, [t, r, [0, 0, 1]]);
}
for (let k = 0; k < 4; k++) { const a = TAU * k / 4 + 0.4; B.cyl([Math.cos(a) * 8, Math.sin(a) * 8, 34], [Math.cos(a), Math.sin(a), 0], 0.5, 16, 6, C.silver); }
for (let k = 0; k < 6; k++) { const a = TAU * k / 6; B.cyl([Math.cos(a) * 4.4, Math.sin(a) * 4.4, 6], [0, 0, 1], 2.4, 24, 10, ACC, C.dark); }
[-1, 1].forEach(s => B.box([s * 19, 0, -26], [34, 0.3, 18], [0.78, 0.79, 0.83, 0.3]));
B.cyl([0, 0, -58], [0, 0, 1], 5, 8, 16, C.dark);
B.dish([0, 0, -70], [0, 0, -1], 7, 10, 24, C.silver, NZ);
});
})();
const SHIP_NOZZLES = { colibri: [[0, 0, -9.6, 1.0]], ranger: [[-1, 0.3, -12.6, 0.6], [0, -1.2, -12.6, 0.6], [1, 0.3, -12.6, 0.6]], atlas: [[0, 1.2, -17, 0.9], [1.4, -0.8, -17, 0.9], [-1.4, -0.8, -17, 0.9]], odisea: [[0, 0, -71, 6.5]] };
const FlightAudio = {
ctx: null, on: true,
init() {
if (this.ctx) { SFX.unlock(); this.active = true; this.setMute(!this.on); return; }
const ctx = SFX.ensure(); if (!ctx) return;
this.ctx = ctx; this.active = true;
const master = this.master = ctx.createGain(); master.gain.value = 0.7; master.connect(SFX.bus.spacecraft);
this.alertOut = ctx.createGain(); this.alertOut.gain.value = 0.7; this.alertOut.connect(SFX.bus.alerts);
const len = ctx.sampleRate * 2, buf = ctx.createBuffer(1, len, ctx.sampleRate), d = buf.getChannelData(0);
let last = 0; for (let i = 0; i < len; i++) { const w = Math.random() * 2 - 1; last = (last + 0.02 * w) / 1.02; d[i] = last * 3.5; }
this.brown = buf;
const wbuf = ctx.createBuffer(1, ctx.sampleRate, ctx.sampleRate), wd = wbuf.getChannelData(0); for (let i = 0; i < wd.length; i++) wd[i] = Math.random() * 2 - 1;
this.white = wbuf;
const src = ctx.createBufferSource(); src.buffer = buf; src.loop = true;
this.engF = ctx.createBiquadFilter(); this.engF.type = 'lowpass'; this.engF.frequency.value = 90;
this.engG = ctx.createGain(); this.engG.gain.value = 0;
src.connect(this.engF); this.engF.connect(this.engG); this.engG.connect(master); src.start();
const hum = ctx.createOscillator(); hum.type = 'sine'; hum.frequency.value = 46; this.humO = hum;
this.humG = ctx.createGain(); this.humG.gain.value = 0; hum.connect(this.humG); this.humG.connect(master); hum.start();
const cab = ctx.createOscillator(); cab.type = 'triangle'; cab.frequency.value = 118;
const cabG = ctx.createGain(); cabG.gain.value = 0.012; cab.connect(cabG); cabG.connect(master); cab.start();
const air = ctx.createBufferSource(); air.buffer = wbuf; air.loop = true;
const airF = ctx.createBiquadFilter(); airF.type = 'bandpass'; airF.frequency.value = 600; airF.Q.value = 0.6;
const airG = ctx.createGain(); airG.gain.value = 0.006; air.connect(airF); airF.connect(airG); airG.connect(master); air.start();
const rcs = ctx.createBufferSource(); rcs.buffer = wbuf; rcs.loop = true;
const rcsF = ctx.createBiquadFilter(); rcsF.type = 'bandpass'; rcsF.frequency.value = 1400; rcsF.Q.value = 1.2;
this.rcsG = ctx.createGain(); this.rcsG.gain.value = 0; rcs.connect(rcsF); rcsF.connect(this.rcsG); this.rcsG.connect(master); rcs.start();
this.setMute(!this.on);
},
setMute(m) {
this.on = !m; if (!this.master) return;
const v = m || !this.active ? 0 : 0.7, t = this.ctx.currentTime;
this.master.gain.setTargetAtTime(v, t, 0.05); this.alertOut.gain.setTargetAtTime(v, t, 0.05);
},
update(thr, rcs) {
if (!this.ctx) return;
const t = this.ctx.currentTime;
this.engG.gain.setTargetAtTime(0.02 + thr * 0.5, t, 0.08);
this.engF.frequency.setTargetAtTime(80 + thr * 380, t, 0.1);
this.humG.gain.setTargetAtTime(thr * 0.12, t, 0.1);
this.humO.frequency.setTargetAtTime(42 + thr * 18, t, 0.2);
this.rcsG.gain.setTargetAtTime(rcs * 0.09, t, 0.03);
},
beep(freq, dur, type, vol) {
if (!this.ctx || !this.on) return;
const t = this.ctx.currentTime, o = this.ctx.createOscillator(), g = this.ctx.createGain();
o.type = type || 'square'; o.frequency.value = freq; g.gain.setValueAtTime(0, t); g.gain.linearRampToValueAtTime(vol || 0.06, t + 0.01); g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
o.connect(g); g.connect(this.alertOut || this.master); o.start(t); o.stop(t + dur + 0.02);
},
ping() { if (!this.ctx) return; const t = this.ctx.currentTime, o = this.ctx.createOscillator(), g = this.ctx.createGain(); o.type = 'sine'; o.frequency.setValueAtTime(660, t); o.frequency.exponentialRampToValueAtTime(990, t + 0.18); g.gain.setValueAtTime(0.0001, t); g.gain.exponentialRampToValueAtTime(0.08, t + 0.02); g.gain.exponentialRampToValueAtTime(0.0001, t + 0.35); o.connect(g); g.connect(this.alertOut || this.master); o.start(t); o.stop(t + 0.4); },
thud() { if (!this.ctx) return; const t = this.ctx.currentTime, s = this.ctx.createBufferSource(), f = this.ctx.createBiquadFilter(), g = this.ctx.createGain(); s.buffer = this.brown; f.type = 'lowpass'; f.frequency.value = 160; g.gain.setValueAtTime(0.8, t); g.gain.exponentialRampToValueAtTime(0.0001, t + 0.9); s.connect(f); f.connect(g); g.connect(this.master); s.start(t); s.stop(t + 1); },
stop() { this.active = false; if (this.master) this.setMute(!this.on); },
};
const Flight = {
on: false, state: null, shipIdx: 0, ship: null, startId: 'tierra', assisted: true,
pos: [0, 0, 0], vel: [0, 0, 0], F: [0, 0, 1], U: [0, 1, 0], w: [0, 0, 0],
stab: true, autoBrake: true, aim: false, match: false, hold: null, orbit: null, intercept: false,
warpIdx: 0, camIdx: 1, free: { az: 0.8, el: 0.25, dist: 3 }, keys: new Set(), stick: null, drag: null,
target: null, ref: null, thr: 0, rcs: 0, accNow: 0, msgs: [], cineT: 0, particles: [], dust: [],
thrCmd: 0, thrAct: 0, thrRef: { v: 0 }, thrAuto: 0, rcsL: 0, rcsU: 0, wsHold: null, cfx: null,
sF: [0, 0, 1], sU: [0, 1, 0], trajT: 0, frameN: 0, alarmT: 0,
L() { return V.cross(this.U, this.F); },
km2s(v) { return V.scale(v, 1000 / AU_KM); },
kmOf(rb) { return V.scale(rb.posS, AU_KM / 1000); },
sRad() { return this.mesh ? this.mesh.R / 1000 : this.ship.len / 2000; },
syncBodies(dtSim) {
for (const rb of World.rb) {
const k = this.kmOf(rb);
if (dtSim > 0 && rb.km) rb.vel = V.scale(V.sub(k, rb.km), 1 / dtSim);
rb.kmPrev = rb.km || k; rb.km = k;
if (!rb.vel) rb.vel = [0, 0, 0];
}
},
initBodies() {
this.grav = World.rb.filter(rb => rb.def.mass && !rb.isCraft && (rb.isSun || rb.def.mass >= 1e21)).map(rb => ({ rb, GM: G_CONST * rb.def.mass * 1e-9 }));
const Msun = World.sun.def.mass;
World.rb.forEach(rb => {
const o = rb.def.orbit; rb.soi = 0;
if (!o || rb.isCraft || !rb.def.mass) return;
const a = o.t === 'jpl' ? o.el[0] * AU_KM : o.t === 'kep' ? o.a * AU_KM : o.a;
const M = rb.parent.isSun ? Msun : rb.parent.def.mass;
if (a && M) rb.soi = Math.max(a * Math.pow(rb.def.mass / M, 0.4), rb.R * 3);
});
World.update(Time.jd - 2 / 86400, 0, 1); this.syncBodies(0);
World.update(Time.jd, 0, 1); this.syncBodies(2);
},
safeDist(rb) { return rb.isCraft ? 1.5 : rb.isSun ? rb.R * 8 : rb.def.mass && rb.R > 100 ? rb.R * 2.2 : Math.max(rb.R * 6, 25); },
refBody() {
let ref = World.sun;
for (const rb of World.rb) if (rb.soi && rb.parent === World.sun && V.dist(this.pos, rb.km) < rb.soi) ref = rb;
if (ref !== World.sun) for (const m of ref.children) if (m.soi && V.dist(this.pos, m.km) < m.soi) ref = m;
return ref;
},
msg(t, kind, ms) { this.msgs.push({ t, kind: kind || 'info', until: performance.now() + (ms || 3500) }); if (kind !== 'quiet') FlightAudio.ping(); },
openHangar() {
if (!this.dom) this.buildDOM();
if (!this.on) {
this.saved = { scale: S.scale, idx: Time.idx, paused: Time.paused, live: Time.live };
UI.stopTour(); if (S.mode !== 'explore') UI.setMode('explore');
UI.select(null); UI.openSection(null); UI.toggleHelp(false);
S.scale = 'real'; ScaleState.set('real', false); World.remapOrbits(); UI.updateScaleBadge();
Time.idx = 0; Time.paused = false;
this.initBodies();
}
this.on = true; this.state = 'hangar';
this.setShip(this.shipIdx); this.placeAtStart();
this.free = { az: 0.9, el: 0.22, dist: 2.2 };
const app = $('#app'); app.classList.add('flying', 'hangar'); app.classList.remove('fl-fly', 'fl-pause');
this.renderHangar();
setTimeout(() => { const b = document.querySelector('#fl-launch'); if (b) b.focus({ preventScroll: true }); }, 60);
PackedModels.preloadShips(this.ship.id);     // modelos del hangar: primero la nave elegida
},
setShip(i) { this.shipIdx = i; this.ship = SHIPS[i]; this.syncModel(); if (this.state === 'hangar') PackedModels.request(this.ship.id, true); },
syncModel() { this.mesh = CraftModels.mesh(this.ship.id) || null; this.nozzles = SHIP_NOZZLES[this.ship.id] || [[0, 0, -this.ship.len / 2, this.ship.len * 0.03]]; },   // provisional mientras se descarga el modelo
placeAtStart() {
const B = World.byId[this.startId], r0 = B.R * (B.id === 'luna' ? 3.5 : 3.2);
const s = V.norm(V.scale(B.km, -1)), dir = V.norm(V.rot(s, [0, 1, 0], 0.65));
this.pos = V.add(B.km, V.scale(dir, r0));
this.startDir = dir; this.startR = r0; this.startB = B;
let F = V.norm(V.rot(V.scale(dir, -1), [0, 1, 0], -0.42));
this.F = F; this.U = V.norm(V.sub([0, 1, 0], V.scale(F, F[1]))); this.w = [0, 0, 0];
this.vel = B.vel.slice();
if (!this.assisted && B.def.mass) {
const t = V.norm(V.cross([0, 1, 0], dir)), vc = Math.sqrt(G_CONST * B.def.mass * 1e-9 / r0);
this.vel = V.add(this.vel, V.scale(t, vc));
}
this.sF = this.F.slice(); this.sU = this.U.slice();
this.ref = this.refBody();
},
launch() {
FlightAudio.init();
this.state = 'fly'; this.placeAtStart();
this.aim = this.match = this.intercept = false; this.hold = null; this.orbit = null; this.warpIdx = 0; this.particles = [];
this.autoBrake = this.assisted; this.stab = true;
this.thrCmd = 0; this.thrAct = 0; this.thrRef = { v: 0 }; this.rcsL = this.rcsU = 0; this.cfx = null;
if (!this.target || this.target === this.startB) this.target = this.startB;
UI.sel = this.target;
const app = $('#app'); app.classList.remove('hangar', 'fl-pause'); app.classList.add('fl-fly');
this.msg('Sistemas en línea. ' + this.ship.name + ' lista para maniobrar.', 'info', 4000);
if (!this.assisted) this.msg('Modo simulación: la gravedad está activa. Partes en una órbita circular.', 'quiet', 6000);
this.refreshBar(); this.updateCockpit();
$('#gl').focus({ preventScroll: true });
},
pause(v) {
if (this.state === 'hangar') return;
this.state = v ? 'pause' : 'fly';
$('#app').classList.toggle('fl-pause', v);
if (v) { FlightAudio.update(0, 0); setTimeout(() => $('#fl-resume').focus({ preventScroll: true }), 30); } else $('#gl').focus({ preventScroll: true });
},
async requestExit() {
if (Settings.state.general.confirmExit && (this.state === 'fly' || this.state === 'pause')) {
const ok = await UI.confirmDialog(I18N.t('¿Salir de la simulación?'), I18N.t('Volverás al explorador y se perderá el vuelo actual.'), I18N.t('Salir'), I18N.t('Seguir volando'));
if (!ok) return;
}
this.exit();
},
exit() {
this.on = false; this.state = null; FlightAudio.update(0, 0); FlightAudio.stop(); Music.setDuck(1);
const app = $('#app'); app.classList.remove('flying', 'hangar', 'fl-fly', 'fl-pause', 'fl-cockpit');
const ref = this.ref || World.byId.tierra;
S.scale = this.saved.scale; ScaleState.set(S.scale, false); World.remapOrbits(); UI.updateScaleBadge();
Time.idx = this.saved.idx; Time.paused = this.saved.paused; Time.live = this.saved.live; if (Time.live) Time.jd = Astro.jdFromDate(new Date()); UI.refreshTime();
World.update(Time.jd, 0, Time.speed);
UI.sel = null; UI.hov = null;
Cam.focus = ref; Cam.follow = true; Cam.target = ref.posS.slice(); Cam.dist = Cam.dDist = Cam.frameDist(ref); Cam.fly = null;
UI.select(ref);
$('#gl').focus({ preventScroll: true });
},
bindInput() {
const cv = $('#gl');
addEventListener('keydown', e => {
if (!this.on) return;
const tag = (e.target.tagName || '').toLowerCase();
if (tag === 'input' || tag === 'select') { if (e.key === 'Escape') e.target.blur(); return; }
if (this.state === 'hangar') { if (e.key === 'Escape') this.exit(); else if (e.key === 'Enter') this.launch(); return; }
if (e.key === 'Escape') { e.preventDefault(); if ($('#fl-targets').classList.contains('open')) this.toggleTargets(false); else this.pause(this.state !== 'pause'); return; }
if (this.state !== 'fly') return;
if ([' ', 'ArrowUp', 'ArrowDown', 'ArrowLeft', 'ArrowRight'].includes(e.key) || e.code === 'Tab') e.preventDefault();
if (!e.repeat) this.onKey(e);
const nk = Keys.norm(e); (this.keyOf = this.keyOf || {})[e.code] = nk; this.keys.add(nk);
});
addEventListener('keyup', e => { this.keys.delete(Keys.norm(e)); if (this.keyOf && this.keyOf[e.code]) this.keys.delete(this.keyOf[e.code]); });
addEventListener('blur', () => this.keys.clear());
cv.addEventListener('pointerdown', e => {
if (!this.on) return;
cv.setPointerCapture(e.pointerId);
this.drag = { x: e.clientX, y: e.clientY, x0: e.clientX, y0: e.clientY, moved: 0, btn: e.button };
if (CAMS[this.camIdx][0] !== 'free' && this.state === 'fly' && e.button === 0) this.stick = { x0: e.clientX, y0: e.clientY, x: 0, y: 0 };
});
cv.addEventListener('pointermove', e => {
if (!this.on || !this.drag) return;
const dx = e.clientX - this.drag.x, dy = e.clientY - this.drag.y;
this.drag.x = e.clientX; this.drag.y = e.clientY; this.drag.moved += Math.abs(dx) + Math.abs(dy);
if (this.stick) { this.stick.x = clamp((e.clientX - this.stick.x0) / 140, -1, 1); this.stick.y = clamp((e.clientY - this.stick.y0) / 140, -1, 1); }
else { this.free.az -= dx * 0.006; this.free.el = clamp(this.free.el + dy * 0.006, -1.5, 1.5); }
});
const up = e => {
if (!this.on || !this.drag) return;
if (this.drag.moved < 5) { const r = cv.getBoundingClientRect(); const rb = World.pick(e.clientX - r.left, e.clientY - r.top); if (rb && this.state === 'fly') this.setTarget(rb); }
this.drag = null; this.stick = null;
};
cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);
cv.addEventListener('wheel', e => { if (!this.on) return; e.preventDefault(); this.free.dist = clamp(this.free.dist * Math.exp(e.deltaY * 0.0012), 0.9, 60); this.zoomK = clamp((this.zoomK || 1) * Math.exp(e.deltaY * 0.0012), 0.5, 8); }, { passive: false });
},
onKey(e) {
const c = e.code;
if (e.shiftKey && /^Digit[1-6]$/.test(c)) return this.setCam(+c.slice(5) - 1);
if (/^Digit[1-4]$/.test(c)) return this.setThrottle([0.25, 0.5, 0.75, 1][+c.slice(5) - 1]);
if (c === 'Digit0') return this.setThrottle(0);
if (c === 'F1' || e.key === '?') { e.preventDefault(); return this.pause(true); }
const key = Keys.norm(e);
switch (Keys.action('flight', key)) {   // acciones reasignables (Ajustes › Controles)
case 'cam': this.cycleCam(1); break;
case 'thr0': this.setThrottle(0); break;
case 'thrUp': this.stepThrottle(1); this.wsHold = { a: 'thrUp', t: 0 }; break;
case 'thrDown': this.stepThrottle(-1); this.wsHold = { a: 'thrDown', t: 0 }; break;
case 'brake': if (this.thrCmd > 0) this.setThrottle(0, 'quiet'); break;
case 'mode': this.toggleAssisted(); break;
case 'stab': this.toggle('stab'); break;
case 'target': e.preventDefault(); this.toggleTargets(); break;
case 'aim': this.toggle('aim'); break;
case 'match': this.toggle('match'); break;
case 'autobrake': this.toggle('autoBrake'); break;
case 'hold': this.toggle('hold'); break;
case 'orbit': this.toggle('orbit'); break;
case 'intercept': this.toggle('intercept'); break;
case 'warpUp': this.setWarp(this.warpIdx + 1); break;
case 'warpDown': this.setWarp(this.warpIdx - 1); break;
case 'sound': this.toggleSound(); break;
}
},
setThrottle(v, src) {
v = clamp(Math.round(v * 100) / 100, DETENTS[0], 1);
if (src !== 'auto' && (this.intercept || this.match || this.hold || this.orbit) && Math.abs(v - this.thrCmd) > 1e-3) {
this.cancelAuto(); this.msg('Control manual del empuje: asistencias de navegación desactivadas.', 'quiet');
}
if (Math.abs(v - this.thrCmd) > 1e-3 && src !== 'quiet') {
FlightAudio.beep(v > this.thrCmd ? 1250 : 980, 0.035, 'triangle', 0.025);
const name = DETENT_NAMES[String(v)];
if (name && src !== 'auto') this.msg('Empuje ' + Math.round(v * 100) + ' %: ' + name.toLowerCase(), 'quiet', 1600);
}
this.thrCmd = v; this.refreshBar();
},
stepThrottle(dir) {
const c = this.thrCmd + 1e-6 * dir;
const next = dir > 0 ? DETENTS.find(d => d > c + 1e-4) : [...DETENTS].reverse().find(d => d < c - 1e-4);
if (next != null) this.setThrottle(next);
},
setCam(i) { this.camIdx = (i + CAMS.length) % CAMS.length; this.cineT = 0; this.updateCockpit(); this.refreshBar(); },
cycleCam(d) { this.setCam(this.camIdx + d); },
updateCockpit() { $('#app').classList.toggle('fl-cockpit', this.state !== 'hangar' && CAMS[this.camIdx][0] === 'cockpit'); },
toggleAssisted() {
this.assisted = !this.assisted; this.autoBrake = this.assisted;
this.msg(this.assisted ? 'Modo asistido: estabilización, freno asistido y sin gravedad.' : 'Modo simulación: inercia y gravedad reales. Mantén la velocidad orbital o caerás.', 'info', 5000);
this.refreshBar();
},
toggleSound() { FlightAudio.setMute(FlightAudio.on); this.refreshBar(); },
toggle(k) {
const T = this.target;
const needT = ['aim', 'match', 'hold', 'intercept'];
if (needT.includes(k) && !T) { this.msg('Elige primero un objetivo (tecla T o clic sobre un objeto).', 'warn'); return; }
if (k === 'stab') { this.stab = !this.stab; this.msg('Estabilización ' + (this.stab ? 'activada' : 'desactivada'), 'quiet'); }
else if (k === 'autoBrake') { this.autoBrake = !this.autoBrake; this.msg('Freno de aproximación ' + (this.autoBrake ? 'activado' : 'desactivado')); }
else if (k === 'aim') { this.aim = !this.aim; if (this.aim) this.msg('Orientando la nave hacia ' + T.def.name); }
else if (k === 'match') { this.cancelAuto('match'); this.match = !this.match; if (this.match) this.msg('Igualando velocidad con ' + T.def.name); if (this.match || this.hold || this.intercept) this.autoEngaged(); }
else if (k === 'hold') { const was = !!this.hold; this.cancelAuto('hold'); this.hold = was ? null : { t: T, d: Math.max(V.dist(this.pos, T.km), this.safeDist(T)) }; if (this.hold) this.msg('Manteniendo distancia con ' + T.def.name); if (this.match || this.hold || this.intercept) this.autoEngaged(); }
else if (k === 'intercept') { const was = this.intercept; this.cancelAuto('intercept'); this.intercept = !was; if (this.intercept) { this.aim = false; this.msg('Piloto automático: rumbo de intercepción a ' + T.def.name); } if (this.match || this.hold || this.intercept) this.autoEngaged(); }
else if (k === 'orbit') {
if (this.orbit) { this.orbit = null; this.msg('Saliendo de órbita'); }
else {
const B = T && T.def.mass && !T.isCraft ? T : this.ref;
const rel = V.sub(this.pos, B.km), r = V.len(rel);
if (B.isSun && !T) { this.msg('Elige un planeta o una luna para orbitar.', 'warn'); return; }
if (r > B.R * 40 && !(B.soi && r < B.soi)) { this.msg('Demasiado lejos de ' + B.def.name + ' para entrar en órbita. Acércate a menos de ' + fmtKm(B.R * 40) + '.', 'warn', 4500); return; }
this.cancelAuto('orbit');
const e1 = V.norm(rel); let n = V.cross(e1, V.sub(this.vel, B.vel)); if (V.len(n) < 1e-9) n = [0, 1, 0]; n = V.norm(n);
const e2 = V.cross(n, e1), GM = G_CONST * B.def.mass * 1e-9, rr = Math.max(r, B.R * 1.15);
this.orbit = { B, r: rr, e1, e2, th: 0, om: Math.sqrt(GM / (rr * rr * rr)) }; this.autoEngaged();
this.msg('Órbita circular asistida alrededor de ' + B.def.name + ': periodo de ' + fmtDur(TAU / this.orbit.om / 86400) + '.', 'info', 5000);
}
}
this.refreshBar();
},
autoEngaged() { if (this.thrCmd !== 0) this.setThrottle(0, 'auto'); },
cancelAuto(keep) { if (keep !== 'match') this.match = false; if (keep !== 'hold') this.hold = null; if (keep !== 'intercept') this.intercept = false; if (keep !== 'orbit') this.orbit = null; },
setTarget(rb) {
if (!rb || rb === this.target) return;
this.target = rb; UI.sel = rb; this.cancelAuto(); this.aim = false;
this.msg('Objetivo: ' + rb.def.name + ' a ' + fmtKm(V.dist(this.pos, rb.km)));
this.refreshBar();
},
setWarp(i) {
const n = clamp(i, 0, WARPS.length - 1);
if (n > this.warpIdx && n > this.maxWarpIdx) { this.msg('Aceleración temporal limitada: hay un cuerpo cerca en tu trayectoria.', 'warn'); return; }
this.warpIdx = n; this.refreshBar();
},
input() {
const k = a => this.keys.has(Keys.code('flight', a)) ? 1 : 0;
const st = this.stick || { x: 0, y: 0 };
return {
l: k('left') - k('right'), u: k('up') - k('down'),
pitch: clamp(k('pitchUp') - k('pitchDown') - st.y, -1, 1), yaw: clamp(k('yawR') - k('yawL') + st.x, -1, 1), roll: k('rollR') - k('rollL'),
brake: k('brake'),
};
},
simulate(dt) {
if (!this.on) return;
if (this.state !== 'fly') {
World.update(Time.jd, dt, 1); this.syncBodies(0);
if (this.state === 'hangar') { this.pos = V.add(this.startB.km, V.scale(this.startDir, this.startR)); this.vel = this.startB.vel.slice(); }
this.ref = this.refBody(); this.thr = 0; this.rcs = 0;
return;
}
const I = this.input(), sh = this.ship;
const manualT = I.l || I.u, manualR = I.pitch || I.yaw || I.roll;
if (manualT && (this.intercept || this.match || this.hold || this.orbit)) { this.cancelAuto(); this.msg('Control manual: asistencias de navegación desactivadas.', 'quiet'); this.refreshBar(); }
if (manualR && this.aim) { this.aim = false; this.refreshBar(); }
if (this.wsHold) {
if (!this.keys.has(Keys.code('flight', this.wsHold.a))) this.wsHold = null;
else { this.wsHold.t += dt; if (this.wsHold.t > 0.4) { this.wsHold.t = 0.22; this.stepThrottle(this.wsHold.a === 'thrUp' ? 1 : -1); } }
}
const sp = SPOOL[sh.id] || [1, 0.8];
const cmd = (this.intercept || this.match || this.hold || this.orbit) ? 0 : this.thrCmd;
this.thrAct = smoothDamp(this.thrAct, cmd, this.thrRef, sp[0], sp[1], dt);
if (Math.abs(this.thrAct) < 1e-4 && cmd === 0) { this.thrAct = 0; this.thrRef.v = 0; }
const autoOn = this.intercept || this.match || this.hold;
this.brk = (this.brk || 0) + ((I.brake ? 1 : 0) - (this.brk || 0)) * (1 - Math.exp(-dt / 0.35));
this.apR = (this.apR || 0) + ((autoOn ? 1 : 0) - (this.apR || 0)) * (1 - Math.exp(-dt / 0.6));
if (this.brk < 0.01 && !I.brake) this.brk = 0;
const rk = 1 - Math.exp(-dt / 0.12);
this.rcsL += (I.l - this.rcsL) * rk; this.rcsU += (I.u - this.rcsU) * rk;
this.updateWarpLimit();
const warp = WARPS[this.warpIdx], dtSim = dt * warp;
Time.idx = Time.speeds.indexOf(warp); Time.paused = false; if (warp > 1) Time.live = false;
Time.jd += dtSim / 86400;
World.update(Time.jd, dt, warp);
this.syncBodies(dtSim);
this.rotate(dt, I);
const vr = V.len(V.sub(this.vel, this.ref.vel));
let hmax = 60;
for (const g of this.grav) { const r = V.dist(this.pos, g.rb.km); if (r < g.rb.R * 60 || g.rb.isSun) hmax = Math.min(hmax, 0.01 * Math.sqrt(r * r * r / g.GM)); }
for (const rb of World.rb) { if (!rb.drawOn && !rb.def.mass) continue; const d = V.dist(this.pos, rb.km) - rb.R; if (d < 1e6) hmax = Math.min(hmax, Math.max(0.02, 0.05 * d / Math.max(vr, 1e-3))); }
let n = Math.min(500, Math.max(1, Math.ceil(dtSim / Math.max(hmax, 1e-3))));
const h = dtSim / n;
this.accNow = 0; this.thr = 0; this.thrAuto = 0; this.rcs = Math.min(1, Math.abs(I.pitch) + Math.abs(I.yaw) + Math.abs(I.roll) + Math.abs(this.rcsL) + Math.abs(this.rcsU));
for (let i = 0; i < n; i++) this.step(h, (i + 1) / n, I);
if (dtSim === 0) this.step(0, 1, I);
this.ref = this.refBody();
this.effects(dt, I);
},
bodyAt(rb, f) { return V.lerp(rb.kmPrev || rb.km, rb.km, f); },
updateWarpLimit() {
let tmin = Infinity;
for (const rb of World.rb) {
if (rb.isCraft && rb !== this.target) continue;
const r = V.sub(rb.km, this.pos), d = V.len(r), clos = V.dot(V.sub(this.vel, rb.vel), r) / Math.max(d, 1e-9);
const free = d - rb.R - this.safeDist(rb) * 0.2;
if (clos > 0) tmin = Math.min(tmin, Math.max(free, 0) / clos);
if (d < rb.R * 3) tmin = Math.min(tmin, 60 * (d / rb.R));
}
let m = WARPS.length - 1;
const k = (this.intercept || this.autoBrake || this.hold || this.match) ? 2.5 : 6;
while (m > 0 && WARPS[m] * k > tmin) m--;
this.maxWarpIdx = m;
if (this.warpIdx > m) { this.warpIdx = m; this.warpLimited = performance.now() + 2500; this.refreshBar(); }
},
rotate(dt, I) {
const sh = this.ship, aA = sh.angAcc * DEG, mA = sh.maxAng * DEG;
const R = V.scale(this.L(), -1), U = this.U, F = this.F;
let des = null, rate = null;
let want = null;
if (this.aim && this.target) want = V.norm(V.sub(this.target.km, this.pos));
if (this.plan && this.plan.dir && (this.intercept || this.match || this.hold)) want = this.plan.dir;
if (want) {
const ax = V.cross(F, want), s = V.len(ax), ang = Math.atan2(s, V.dot(F, want));
const axis = s > 1e-9 ? V.scale(ax, 1 / s) : U;
const sp = Math.min(mA, Math.sqrt(2 * aA * ang * 0.8));
const wv = V.scale(axis, sp);
des = [V.dot(wv, R), V.dot(wv, U), 0];
}
const inp = [I.pitch, -I.yaw, I.roll];
for (let k = 0; k < 3; k++) {
let target;
if (inp[k]) target = this.assisted ? inp[k] * mA : null;
else if (des) target = des[k];
else if (this.stab || this.assisted) target = 0;
else target = null;
if (target == null) this.w[k] = clamp(this.w[k] + inp[k] * aA * dt, -mA * 1.5, mA * 1.5);
else this.w[k] += clamp(target - this.w[k], -aA * dt, aA * dt);
}
const wv = V.add(V.add(V.scale(R, this.w[0]), V.scale(U, this.w[1])), V.scale(F, this.w[2]));
const om = V.len(wv);
if (om > 1e-12) { const ax = V.scale(wv, 1 / om); this.F = V.rot(this.F, ax, om * dt); this.U = V.rot(this.U, ax, om * dt); }
this.F = V.norm(this.F); this.U = V.norm(V.sub(this.U, V.scale(this.F, V.dot(this.U, this.F))));
},
thrustToward(D, h, power) {
const sh = this.ship, main = sh.accel / 1000 * (power || 1), m = V.len(D);
if (m < 1e-9) return [0, 0, 0];
const d = V.scale(D, 1 / m), al = V.dot(this.F, d);
let a = al > 0.985 ? main : al < -0.985 ? main * 0.4 : main * 0.25;
if (this.assisted) a = Math.max(a, main * 0.6);
a = Math.min(a, h > 0 ? m / h : a);
return V.scale(d, a);
},
step(h, f, I) {
const sh = this.ship, main = sh.accel / 1000, L = this.L();
let acc = [0, 0, 0], thr = 0;
const T = this.target;
if (this.thrAct !== 0) { acc = V.add(acc, V.scale(this.F, this.thrAct * main)); thr = Math.max(0, this.thrAct); }
if (Math.abs(this.rcsL) > 1e-3 || Math.abs(this.rcsU) > 1e-3) acc = V.add(acc, V.add(V.scale(L, this.rcsL * main * 0.25), V.scale(this.U, this.rcsU * main * 0.25)));
const ref = this.ref, vRef = ref.vel, vrel = V.sub(this.vel, vRef);
this.plan = null;
if (this.orbit) {
const O = this.orbit, B = this.bodyAt(O.B, f);
O.th += O.om * h;
const c = Math.cos(O.th), s = Math.sin(O.th);
this.pos = V.add(B, V.add(V.scale(O.e1, O.r * c), V.scale(O.e2, O.r * s)));
this.vel = V.add(O.B.vel, V.scale(V.add(V.scale(O.e1, -s), V.scale(O.e2, c)), O.r * O.om));
this.plan = { dir: null };
return;
}
if (this.brk > 0.01) {
const b0 = this.assisted ? V.scale(V.norm(vrel), -Math.min(main, h > 0 ? V.len(vrel) / h : main)) : this.thrustToward(V.scale(vrel, -1), h);
const b = V.scale(b0, this.brk);
if (V.len(vrel) > 1e-6) { acc = V.add(acc, b); thr = Math.max(thr, 0.6 * this.brk); }
}
if (T && (this.intercept || this.match || this.hold)) {
const Tp = this.bodyAt(T, f), r = V.sub(Tp, this.pos), d = V.len(r), rh = V.scale(r, 1 / Math.max(d, 1e-9));
const vr = V.sub(this.vel, T.vel), aPlan = sh.accel / 1000 * 0.35;
let vdes = [0, 0, 0];
if (this.intercept) {
const ds = this.safeDist(T), sp = Math.min(sh.maxSpeed * 0.97, Math.sqrt(2 * aPlan * Math.max(d - ds, 0)));
vdes = V.scale(rh, sp);
if (d < ds * 1.3 && V.len(vr) < Math.max(0.005, sp * 0.2 + 0.002)) { this.intercept = false; this.hold = { t: T, d: Math.max(d, ds) }; this.arrived = T; }
} else if (this.hold) vdes = V.scale(rh, clamp((d - this.hold.d) * 0.05, -sh.maxSpeed * 0.2, sh.maxSpeed * 0.2));
const D = V.sub(vdes, vr);
const a = V.scale(this.thrustToward(D, h), Math.max(this.apR || 0, 0.05));
acc = V.add(acc, a);
if (V.len(a) > 1e-9) { const al = V.dot(a, this.F) / main; this.thrAuto = Math.max(this.thrAuto, clamp(al, 0, 1)); thr = Math.max(thr, clamp(al, 0, 1)); }
this.plan = { dir: V.len(D) > sh.accel / 1000 * 2 ? V.norm(D) : rh };
}
if (this.assisted && Math.abs(this.rcsL) < 0.05 && Math.abs(this.rcsU) < 0.05 && !this.intercept && !this.match && !this.hold && !I.brake) {
const lat = V.sub(vrel, V.scale(this.F, V.dot(vrel, this.F)));
const lm = V.len(lat);
if (lm > 1e-7) acc = V.add(acc, V.scale(lat, -Math.min(main * 0.25, h > 0 ? lm / h : 0) / lm));
}
this.braking = false;
if (this.autoBrake) {
for (const rb of [this.target, ref].filter(Boolean)) {
const r = V.sub(this.bodyAt(rb, f), this.pos), d = V.len(r), rh = V.scale(r, 1 / Math.max(d, 1e-9));
const clos = V.dot(V.sub(this.vel, rb.vel), rh), free = d - rb.R - this.safeDist(rb) * 0.5;
const a = main * 0.9;
if (clos > 0 && clos * clos > 2 * a * 0.8 * Math.max(free, 0)) { acc = V.add(acc, V.scale(rh, -Math.min(a, h > 0 ? clos / h : a))); this.braking = true; thr = Math.max(thr, 0.5); }
}
}
let g = [0, 0, 0];
if (!this.assisted) for (const G of this.grav) {
const r = V.sub(this.bodyAt(G.rb, f), this.pos), d2 = V.dot(r, r), d = Math.sqrt(d2);
if (d < G.rb.R * 0.5) continue;
g = V.add(g, V.scale(r, G.GM / (d2 * d)));
}
const tA = V.len(acc);
this.accNow = Math.max(this.accNow, tA * 1000);
this.thr = Math.max(this.thr, Math.min(1, thr));
this.vel = V.add(this.vel, V.scale(V.add(acc, g), h));
const vl = V.sub(this.vel, vRef), sp = V.len(vl);
if (sp > sh.maxSpeed) this.vel = V.add(vRef, V.scale(vl, sh.maxSpeed / sp));
this.pos = V.add(this.pos, V.scale(this.vel, h));
for (const rb of World.rb) {
if (rb.isCraft) continue;
const B = this.bodyAt(rb, f), r = V.sub(this.pos, B), d = V.len(r), lim = rb.R + this.sRad();
if (d < lim && rb.R > 0) {
this.pos = V.add(B, V.scale(r, (lim * 1.0005) / d)); this.vel = rb.vel.slice(); this.cancelAuto(); this.warpIdx = 0; this.thrCmd = 0; this.thrAct = 0; this.thrRef.v = 0;
this.msg('Contacto con ' + rb.def.name + '. La nave se ha detenido sobre la superficie.', 'alert', 5000); FlightAudio.thud(); this.refreshBar();
}
}
},
effects(dt, I) {
const sh = this.ship, len = sh.len / 1000;
if (this.thr > 0.02) for (const nz of this.nozzles) for (let k = 0; k < 3; k++) {
if (Math.random() > this.thr) continue;
const j = [(Math.random() - 0.5) * nz[3] * 0.6, (Math.random() - 0.5) * nz[3] * 0.6, 0];
this.particles.push({ p: [nz[0] + j[0], nz[1] + j[1], nz[2]], v: [j[0] * 1.5, j[1] * 1.5, -(sh.len * 2.2) * (0.7 + Math.random() * 0.6)], age: 0, life: 0.5 + Math.random() * 0.4 });
}
this.particles = this.particles.filter(p => (p.age += dt) < p.life);
this.particles.forEach(p => { p.p = V.add(p.p, V.scale(p.v, dt)); });
if (this.particles.length > 260) this.particles.splice(0, this.particles.length - 260);
const now = performance.now();
this.alert = null;
for (const rb of [this.target, this.ref].filter(Boolean)) {
const r = V.sub(rb.km, this.pos), d = V.len(r), clos = V.dot(V.sub(this.vel, rb.vel), r) / Math.max(d, 1e-9);
const tti = clos > 0 ? (d - rb.R) / clos : Infinity;
if (tti < 20 * Math.max(1, WARPS[this.warpIdx] / 1)) { this.alert = { kind: 'prox', txt: 'Proximidad: impacto con ' + rb.def.name + ' en ' + fmtDur(tti / 86400 / WARPS[this.warpIdx] * WARPS[this.warpIdx]) }; break; }
if (d - rb.R < this.safeDist(rb) && clos > Math.sqrt(2 * sh.accel / 1000 * 0.5 * Math.max(d - rb.R, 0.001))) { this.alert = { kind: 'speed', txt: 'Velocidad excesiva para la aproximación a ' + rb.def.name }; break; }
}
if (this.alert && now > this.alarmT) { this.alarmT = now + (this.alert.kind === 'prox' ? 600 : 1100); FlightAudio.beep(this.alert.kind === 'prox' ? 880 : 620, 0.16, 'square', 0.05); if (this.alert.kind === 'speed') setTimeout(() => FlightAudio.beep(480, 0.16, 'square', 0.05), 180); }
if (this.arrived) { this.msg('Has llegado a ' + this.arrived.def.name + '. Manteniendo distancia.', 'info', 5000); this.arrived = null; this.refreshBar(); }
FlightAudio.update(this.state === 'fly' ? this.thr : 0, this.state === 'fly' ? Math.min(1, this.rcs + (this.braking ? 0.3 : 0)) : 0);
Music.setDuck(this.state === 'fly' && FlightAudio.on ? 0.55 : 1);
},
camera(dt) {
const sh = this.ship, len = sh.len / 1000, F = this.F, U = this.U, L = this.L(), p = this.pos;
const mode = this.state === 'hangar' ? 'free' : CAMS[this.camIdx][0], C = CAMFX[mode];
const k = 1 - Math.exp(-dt * 3.2);
this.sF = V.norm(V.lerp(this.sF, F, k)); this.sU = V.norm(V.lerp(this.sU, V.sub(U, V.scale(this.sF, V.dot(U, this.sF))), k));
const fx = this.cfx || (this.cfx = { acc: 0, accR: { v: 0 }, prev: 0, env: 0, back: 0, backR: { v: 0 }, fov: C.fov, fovR: { v: 0 }, t: 0, ph: [0, 1, 2, 3, 4, 5].map(() => Math.random() * TAU), fz: { az: this.free.az, el: this.free.el, d: this.free.dist }, pw: 0 });
const accN = this.state === 'fly' ? clamp(this.accNow / Math.max(sh.accel, 1e-6), 0, 1) : 0;
fx.acc = smoothDamp(fx.acc, accN, fx.accR, 0.45, 4, dt);
const jerk = dt > 0 ? Math.abs(fx.acc - fx.prev) / dt : 0; fx.prev = fx.acc;
const ign = this.thrAct > 0.03 && fx.pw <= 0.03 ? 0.35 : 0; fx.pw = this.thrAct;
const man = clamp((Math.abs(this.rcsL) + Math.abs(this.rcsU)) * 0.25 + Math.hypot(...this.w) / Math.max(sh.maxAng * DEG, 1e-6) * 0.08, 0, 0.3);
fx.env = Math.max(fx.env * Math.exp(-dt / 0.5), clamp(jerk * 0.55, 0, 1), ign, man * 0.5);
const backT = fx.acc * C.back * (UI.reducedMotion ? 0 : 1);
fx.back = smoothDamp(fx.back, backT, fx.backR, backT > fx.back ? 0.9 : 1.8, 2, dt);
const still = UI.reducedMotion ? 0 : 1;
const fovT = C.fov + fx.acc * C.fovUp * still;
fx.fov = smoothDamp(fx.fov, fovT, fx.fovR, 0.8, 20, dt);
fx.t += dt;
const amp = (fx.env * C.shake + Math.max(this.thrAct, 0) * C.steady) * still;
const nz = i => Math.sin(fx.t * 11.3 + fx.ph[i]) * 0.55 + Math.sin(fx.t * 19.7 + fx.ph[i + 3]) * 0.3 + Math.sin(fx.t * 31.1 + fx.ph[i] * 2.1) * 0.15;
let pos, look, up = U;
const z = this.zoomK || 1, bk = 1 + fx.back;
if (mode === 'cockpit') {
pos = V.add(p, V.add(V.scale(F, len * (0.36 - C.back * fx.acc)), V.scale(U, len * 0.07)));
const yaw = nz(0) * amp * 0.6, pit = nz(1) * amp * 0.6;
look = V.add(pos, V.norm(V.add(F, V.add(V.scale(L, yaw), V.scale(U, pit)))));
} else if (mode === 'chase') {
const sF = this.sF, sU = this.sU;
pos = V.add(p, V.add(V.scale(sF, -len * 2.3 * z * bk), V.scale(sU, len * 0.55 * z * (1 + fx.back * 0.3))));
look = V.add(p, V.scale(sF, len * 0.9));
up = sU;
} else if (mode === 'front') { pos = V.add(p, V.add(V.scale(F, len * 2.6 * z * bk), V.scale(U, len * 0.35 * z))); look = p; }
else if (mode === 'side') { pos = V.add(p, V.add(V.scale(L, len * 2.8 * z * bk), V.scale(U, len * 0.3 * z))); look = p; }
else if (mode === 'cine') {
this.cineT += dt; const shot = Math.floor(this.cineT / 8) % 4, s = easeInOut((this.cineT % 8) / 8);
const offs = [[-3.2 + s, 1.4 - s * 2.2, 0.3], [3.6 - s * 0.8, -2.2, 0.8 + s * 0.5], [-0.6, 4.5 - s * 1.5, 1.2], [-6 - s * 2, 0.5, 2.6 - s]];
const o = offs[shot], sF = this.sF, sU = this.sU, sL = V.cross(sU, sF);
pos = V.add(p, V.scale(V.add(V.add(V.scale(sF, o[0] - fx.back * 2), V.scale(sL, o[1])), V.scale(sU, o[2])), len * z));
look = V.add(p, V.scale(sF, len * 0.4)); up = [0, 1, 0];
} else {
if (this.state === 'hangar') this.free.az += dt * 0.12;
const kf = 1 - Math.exp(-dt * 9);
fx.fz.az += wrapPi(this.free.az - fx.fz.az) * kf; fx.fz.el += (this.free.el - fx.fz.el) * kf; fx.fz.d += (this.free.dist - fx.fz.d) * kf;
const ce = Math.cos(fx.fz.el), d = fx.fz.d * len * bk;
pos = V.add(p, V.scale([ce * Math.sin(fx.fz.az), Math.sin(fx.fz.el), ce * Math.cos(fx.fz.az)], d)); look = p; up = [0, 1, 0];
}
if (mode !== 'cockpit') {
const f0 = V.norm(V.sub(look, pos)), r0 = V.norm(V.cross(f0, Math.abs(V.dot(f0, up)) > 0.995 ? [1, 0, 0] : up)), u0 = V.cross(r0, f0);
const off = V.add(V.scale(r0, nz(0) * amp * len), V.scale(u0, nz(1) * amp * len));
pos = V.add(pos, off); look = V.add(look, V.scale(off, 0.5));
}
const camS = this.km2s(pos), f = V.norm(V.sub(look, pos));
World.cam = camS; Cam.pos = camS;
World.view = M4.viewRot(f, Math.abs(V.dot(f, up)) > 0.995 ? [1, 0, 0] : up, World.view);
const fov = this.state === 'hangar' ? 50 : fx.fov;
World.proj = M4.perspInf(fov * DEG, GLX.W / GLX.H, Math.max(len * 0.004 * 1000 / AU_KM, 1e-12), World.proj);
if (this.state === 'hangar' && innerWidth >= 760) World.proj[8] = (396 - 316) / innerWidth;
if (this.state === 'hangar' && innerWidth < 760) World.proj[9] = -0.25;
this.camMode = mode; this.camKm = pos; this.fovNow = fov;
},
frame(dt, t) { this.simulate(dt); this.renderFrame(t, dt); },
renderFrame(t, dt) {
this.camera(dt || 0.016);
App.prepare();
World.render(t); World.frame++;
UI.updateLabels();
this.hud();
},
drawShip(useP) {
if (this.camMode === 'cockpit') return;
this.syncModel();
if (!this.mesh) return;                       // el modelo aún se está descargando
const tex = !!this.mesh.groups, P = tex ? World.P.craftTex : World.P.craft, cam = World.cam;
const rel = V.sub(this.km2s(this.pos), cam), s = this.mesh.R / AU_KM;
const m = M4.fromBasis(this.L(), this.U, this.F, s, rel);
const pr = useP(P), ob = new Float32Array(16); let n = 0;
for (const rb of [this.ref, World.byId.tierra, World.byId.luna]) { if (!rb || rb.isSun || n >= 4) continue; const r2 = V.sub(rb.posS, cam); ob[n * 4] = r2[0]; ob[n * 4 + 1] = r2[1]; ob[n * 4 + 2] = r2[2]; ob[n * 4 + 3] = rb.rS; n++; }
const nz = this.nozzles[0], ep = this.km2s(V.add(this.pos, this.bodyVec([nz[0], nz[1], nz[2] - this.ship.len * 0.12])));
GLX.setAll(pr, { u_model: m, u_occ: ob, u_occN: n, u_sunR: World.shadowSunR, u_ambient: S.ambient + 0.01, u_hover: 0, u_eng: [...V.sub(ep, cam), this.thr * 1.3], u_engCol: this.ship.eng, u_engR: this.ship.len * 0.35 / AU_KM });
if (tex) World.drawGroups(this.mesh, pr); else GLX.draw(this.mesh, pr);
},
bodyVec(v) { const L = this.L(); return V.scale(V.add(V.add(V.scale(L, v[0]), V.scale(this.U, v[1])), V.scale(this.F, v[2])), 1 / 1000); },
drawFx(useP) {
const gl = GLX.gl, cam = World.cam, sh = this.ship, P = World.P, lenS = sh.len / AU_KM;
if (!this.fx) {
this.fx = {
pts: GLX.mesh({ a_pos: { data: new Float32Array(600 * 3), size: 3 }, a_col: { data: new Float32Array(600 * 4), size: 4 } }, null, gl.POINTS, true),
traj: GLX.mesh({ a_pos: { data: new Float32Array(241 * 3), size: 3 }, a_frac: { data: new Float32Array(241).map((_, i) => 1 - i / 241), size: 1 } }, null, gl.LINE_STRIP, true),
pb: new Float32Array(600 * 3), cb: new Float32Array(600 * 4),
};
this.dust = Array.from({ length: 320 }, () => [Math.random(), Math.random(), Math.random()]);
}
if (this.state === 'hangar') return;
this.trajT -= 1;
if (this.trajT <= 0 || !this.trajPts) { this.predict(); this.trajT = 6; }
if (this.trajPts) {
const b = new Float32Array(241 * 3), refS = this.km2s(this.ref.km);
this.trajPts.forEach((q, i) => { const s = V.sub(V.add(refS, this.km2s(q)), cam); b[i * 3] = s[0]; b[i * 3 + 1] = s[1]; b[i * 3 + 2] = s[2]; });
GLX.update(this.fx.traj, 'a_pos', b);
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'alpha' });
const pr = useP(P.orbit);
GLX.setAll(pr, { u_model: M4.ident(), u_color: lin('#86d3c6'), u_alpha: 0.75, u_cur: 1, u_fade: 0.85 });
GLX.draw(this.fx.traj, pr, this.trajPts.length);
}
GLX.state({ depthTest: true, depthWrite: false, cull: 'none', blend: 'add' });
const pg = useP(P.glow), thr = this.thr;
for (const nz of this.nozzles) {
const c0 = V.sub(this.km2s(V.add(this.pos, this.bodyVec([nz[0], nz[1], nz[2]]))), cam);
const r = nz[3] / AU_KM;
GLX.setAll(pg, { u_center: c0, u_size: r * (1.8 + thr * 1.4), u_color: sh.eng, u_int: (0.08 + thr * 0.7) * Gfx.flash(), u_sharp: 7, u_near: 0 });
GLX.draw(World.meshes.quad, pg);
if (thr > 0.05) for (let k = 1; k <= 5; k++) {
const c = V.sub(this.km2s(V.add(this.pos, this.bodyVec([nz[0], nz[1], nz[2] - nz[3] * 2.2 * k * (0.6 + thr)]))), cam);
GLX.setAll(pg, { u_center: c, u_size: r * (1.3 + k * 0.45) * (0.6 + thr * 0.7), u_color: sh.eng, u_int: thr * 0.38 / k, u_sharp: 6, u_near: 0 });
GLX.draw(World.meshes.quad, pg);
}
}
if (this.rcs > 0.05 || this.braking) {
const pts = [[1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1], [0, 1, 0.8], [0, -1, -0.8]];
pts.forEach(q => {
if (Math.random() > 0.65) return;
const c = V.sub(this.km2s(V.add(this.pos, this.bodyVec([q[0] * sh.len * 0.12, q[1] * sh.len * 0.06, q[2] * sh.len * 0.35]))), cam);
GLX.setAll(pg, { u_center: c, u_size: lenS * 0.05, u_color: [0.9, 0.95, 1.0], u_int: 0.5, u_sharp: 9, u_near: 0 });
GLX.draw(World.meshes.quad, pg);
});
}
const pb = this.fx.pb, cb = this.fx.cb; let n = 0;
for (const p of this.particles) {
if (n >= 600) break;
const s = V.sub(this.km2s(V.add(this.pos, this.bodyVec(p.p))), cam), a = (1 - p.age / p.life);
pb.set(s, n * 3); const e = sh.eng; cb.set([e[0] * a * 1.4, e[1] * a * 1.4, e[2] * a * 1.4, a], n * 4); n++;
}
const ref = this.ref, vrel = V.len(V.sub(this.vel, ref.vel));
const target = clamp(vrel * 1.6, 0.06, 4e5);
this.dustS = this.dustS ? this.dustS * Math.exp(-0.05) + target * (1 - Math.exp(-0.05)) : target;
const Sz = this.dustS, c = V.sub(this.camKm, ref.km);
const dustN = Settings.state.graphics.dust ? Math.round(this.dust.length * Gfx.particles()) : 0;
for (let di = 0; di < dustN; di++) {
const d = this.dust[di];
if (n >= 600) break;
const q = [0, 1, 2].map(i => { const o = d[i] * Sz; return o + Sz * Math.round((c[i] - o) / Sz); });
const rel = V.sub(q, c), dd = V.len(rel) / (Sz * 0.5);
if (dd > 1) continue;
const a = (1 - dd) * 0.35 * clamp(vrel / 0.2, 0.15, 1);
pb.set(this.km2s(rel), n * 3); cb.set([0.7 * a, 0.78 * a, 0.9 * a, a], n * 4); n++;
}
if (n) {
GLX.update(this.fx.pts, 'a_pos', pb); GLX.update(this.fx.pts, 'a_col', cb);
const pm = useP(P.mark); GLX.setAll(pm, { u_dpr: GLX.dpr, u_size: 3 });
GLX.draw(this.fx.pts, pm, n);
}
},
predict() {
const ref = this.ref, GMr = ref.def.mass ? G_CONST * ref.def.mass * 1e-9 : 0;
let p = V.sub(this.pos, ref.km), v = V.sub(this.vel, ref.vel);
const r0 = V.len(p), sp = Math.max(V.len(v), 1e-6);
let T;
if (!this.assisted && GMr && sp * sp / 2 - GMr / r0 < 0) { const a = 1 / (2 / r0 - sp * sp / GMr); T = TAU * Math.sqrt(a * a * a / GMr) * 1.02; }
else T = this.target ? clamp(1.3 * V.dist(this.pos, this.target.km) / sp, 60, 3e7) : clamp(r0 * 3 / sp, 60, 3e7);
T = Math.min(T, 3e7);
const pts = [p.slice()], N = 240, sub = 6, h = T / N / sub;
for (let i = 0; i < N; i++) {
for (let k = 0; k < sub; k++) {
if (!this.assisted && GMr) { const d2 = V.dot(p, p), d = Math.sqrt(d2); if (d < ref.R) break; v = V.add(v, V.scale(p, -GMr / (d2 * d) * h)); }
p = V.add(p, V.scale(v, h));
}
pts.push(p.slice());
}
this.trajPts = pts;
},
buildDOM() {
const el = $('#flight');
el.innerHTML = `
<svg id="cockpit" viewBox="0 0 1600 900" preserveAspectRatio="none" aria-hidden="true">
<path fill-rule="evenodd" fill="#04060c" d="M0 0H1600V900H0Z M120 70 L1480 70 L1560 330 L1500 760 L1080 840 L520 840 L100 760 L40 330 Z"/>
<path fill="none" stroke="rgba(143,195,255,.18)" stroke-width="2" d="M120 70 L1480 70 L1560 330 L1500 760 L1080 840 L520 840 L100 760 L40 330 Z M800 70 L800 40 M520 840 L560 900 M1080 840 L1040 900"/>
<path fill="#04060c" d="M770 70 L830 70 L812 120 L788 120 Z"/>
</svg>
<div id="fl-mk" aria-hidden="true">
<div id="mk-bore" class="mk"></div><div id="mk-pro" class="mk"></div><div id="mk-retro" class="mk"></div>
<div id="mk-tgt" class="mk"><i></i><div class="t"><b></b><span></span></div></div><div id="mk-arrow" class="mk"></div>
</div>
<div id="hud-top" class="hud"><span id="h-ship"></span><span id="h-mode"></span><span id="h-cam"></span><span id="h-warp"></span></div>
<div id="hud-warn" class="hud" role="status" aria-live="polite"></div>
<section id="hud-l" class="hud hud-col" aria-label="Navegación">
<h4>Navegación</h4><dl id="h-nav"></dl>
</section>
<section id="hud-thr" class="hud" aria-label="Acelerador">
<div class="thr-head"><span>Empuje</span><b id="thr-act">0 %</b></div>
<div class="thr-body">
<div id="thr-track" role="slider" tabindex="-1" aria-label="Empuje solicitado" aria-valuemin="-25" aria-valuemax="100" aria-valuenow="0">
<span class="thr-tick" style="bottom:0.00%"><em>−25</em></span><span class="thr-tick z" style="bottom:20.00%"><em>0</em></span><span class="thr-tick" style="bottom:28.00%"><em>10</em></span><span class="thr-tick" style="bottom:40.00%"><em>25</em></span><span class="thr-tick" style="bottom:60.00%"><em>50</em></span><span class="thr-tick" style="bottom:80.00%"><em>75</em></span><span class="thr-tick" style="bottom:92.00%"><em>90</em></span><span class="thr-tick" style="bottom:100.00%"><em>100</em></span><i id="thr-fill"></i><i id="thr-lever"></i>
</div>
<div class="thr-btns">
<button data-thr="up" aria-label="Aumentar empuje" data-tip="Aumentar empuje (W)">+</button>
<button data-thr="1" data-tip="Empuje 100 % (4)">100</button><button data-thr="0.75" data-tip="Empuje 75 % (3)">75</button><button data-thr="0.5" data-tip="Empuje 50 % (2)">50</button><button data-thr="0.25" data-tip="Empuje 25 % (1)">25</button><button data-thr="0" data-tip="Motores en reposo (X)">0</button>
<button data-thr="down" aria-label="Disminuir empuje" data-tip="Disminuir empuje (S)">−</button>
</div>
</div>
<div class="thr-foot"><span>Solicitado <b id="thr-cmd">0 %</b></span><span id="thr-acc">0 m/s²</span><span id="thr-name"></span></div>
</section>
<section id="hud-r" class="hud hud-col" aria-label="Objetivo">
<h4>Objetivo</h4><dl id="h-tgt"></dl>
</section>
<div id="hud-scale" class="hud">Escala real: distancias, tamaños, velocidades y tiempos coherentes. Las naves son conceptuales y se ignora la relatividad.</div>
<div id="fl-bar" role="toolbar" aria-label="Controles de vuelo"></div>
<div id="fl-targets" role="dialog" aria-label="Elegir objetivo"><div class="ft-head"><input id="ft-q" type="search" placeholder="Buscar objetivo…" aria-label="Buscar objetivo"><button class="txt-btn" id="ft-close">Cerrar</button></div><div id="ft-list"></div></div>
<div id="fl-pause" role="dialog" aria-modal="true" aria-labelledby="fp-title"><div class="fp-card">
<h2 id="fp-title">Vuelo en pausa</h2>
<div class="fp-grid">
<dl><dt><kbd>W</kbd> <kbd>S</kbd></dt><dd>Subir y bajar el empuje por posiciones</dd><dt><kbd>1</kbd>–<kbd>4</kbd></dt><dd>Empuje 25, 50, 75 y 100 %</dd><dt><kbd>X</kbd> o <kbd>0</kbd></dt><dd>Motores en reposo</dd><dt><kbd>A</kbd> <kbd>D</kbd></dt><dd>Propulsores laterales</dd><dt><kbd>R</kbd> <kbd>F</kbd></dt><dd>Propulsores arriba y abajo</dd><dt><kbd>↑</kbd> <kbd>↓</kbd> o arrastrar</dt><dd>Cabeceo</dd><dt><kbd>←</kbd> <kbd>→</kbd> o arrastrar</dt><dd>Guiñada</dd><dt><kbd>Q</kbd> <kbd>E</kbd></dt><dd>Alabeo</dd><dt><kbd>Espacio</kbd></dt><dd>Frenar (corta el empuje)</dd></dl>
<dl><dt><kbd>Z</kbd></dt><dd>Estabilización automática</dd><dt><kbd>V</kbd></dt><dd>Modo asistido o simulación</dd><dt><kbd>C</kbd> o <kbd>Mayús</kbd> + <kbd>1</kbd>–<kbd>6</kbd></dt><dd>Cámara</dd><dt><kbd>,</kbd> <kbd>.</kbd></dt><dd>Aceleración temporal</dd><dt><kbd>T</kbd> o clic</dt><dd>Elegir objetivo</dd><dt><kbd>M</kbd></dt><dd>Silenciar sonido</dd><dt><kbd>Esc</kbd></dt><dd>Pausa</dd></dl>
<dl><dt><kbd>G</kbd></dt><dd>Orientar hacia el objetivo</dd><dt><kbd>B</kbd></dt><dd>Igualar velocidad</dd><dt><kbd>N</kbd></dt><dd>Freno de aproximación</dd><dt><kbd>H</kbd></dt><dd>Mantener distancia</dd><dt><kbd>O</kbd></dt><dd>Entrar o salir de órbita</dd><dt><kbd>I</kbd></dt><dd>Interceptar (piloto automático)</dd></dl>
</div>
<p class="note">En modo asistido no hay gravedad y la nave corrige su deriva. En modo simulación la nave conserva su inercia y la gravedad del Sol, los planetas y las lunas grandes curva tu trayectoria.</p>
<div class="fp-actions"><button class="cta" id="fl-resume">Continuar vuelo</button><button class="txt-btn" id="fl-tohangar">Cambiar de nave</button><button class="txt-btn" id="fl-exit">Salir del modo de vuelo</button></div>
</div></div>
<div id="hangar" role="dialog" aria-labelledby="hg-title">
<div class="hg-panel">
<button class="hg-home" id="hg-home"><svg viewBox="0 0 24 24" width="16" height="16" aria-hidden="true"><path d="M15 5l-7 7 7 7" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>Volver al inicio</button>
<h2 id="hg-title">Elige tu nave</h2>
<p class="note">Naves conceptuales: sus prestaciones son teóricas o ficticias. El Sistema Solar, en cambio, está a escala real.</p>
<div id="hg-list" role="radiogroup" aria-label="Naves"></div>
</div>
<div class="hg-info">
<span id="hg-type"></span><h3 id="hg-name"></h3><p id="hg-desc"></p>
<dl id="hg-stats"></dl>
<div id="hg-extra"></div>
<div class="hg-opts">
<label>Punto de partida <select id="hg-start">${STARTS.map(([id, t]) => `<option value="${id}">${t}</option>`).join('')}</select></label>
<div class="seg" role="radiogroup" aria-label="Modo de vuelo"><button role="radio" data-as="1">Asistido</button><button role="radio" data-as="0">Simulación</button></div>
</div>
<p id="hg-credit" class="note"></p>
<div class="hg-actions"><button class="cta" id="fl-launch">Despegar</button><button class="txt-btn" id="hg-back">Volver al explorador</button></div>
<div id="hg-brand" hidden><img alt="" decoding="async"></div>
</div>
</div>`;
this.dom = true;
$('#fl-resume').addEventListener('click', () => this.pause(false));
$('#fl-exit').addEventListener('click', () => this.requestExit());
$('#fl-tohangar').addEventListener('click', () => { this.pause(false); this.openHangar(); });
$('#fl-launch').addEventListener('click', () => this.launch());
$('#hg-back').addEventListener('click', () => this.exit());
$('#hg-home').addEventListener('click', () => App.home());
$('#hg-start').addEventListener('change', e => { this.startId = e.target.value; this.target = World.byId[this.startId]; this.placeAtStart(); });
document.querySelectorAll('[data-as]').forEach(b => b.addEventListener('click', () => { this.assisted = b.dataset.as === '1'; this.renderHangar(); }));
$('#ft-close').addEventListener('click', () => this.toggleTargets(false));
$('#ft-q').addEventListener('input', () => this.renderTargets());
document.querySelectorAll('[data-thr]').forEach(b => b.addEventListener('click', () => {
const a = b.dataset.thr; if (a === 'up') this.stepThrottle(1); else if (a === 'down') this.stepThrottle(-1); else this.setThrottle(+a);
$('#gl').focus({ preventScroll: true });
}));
const tr = $('#thr-track');
const fromY = e => { const r = tr.getBoundingClientRect(); let v = (1 - (e.clientY - r.top) / r.height) * 1.25 - 0.25; const near = DETENTS.find(d => Math.abs(d - v) < 0.035); return clamp(near != null ? near : v, -0.25, 1); };
tr.addEventListener('pointerdown', e => { tr.setPointerCapture(e.pointerId); this.thrDrag = true; this.setThrottle(fromY(e)); });
tr.addEventListener('pointermove', e => { if (this.thrDrag) this.setThrottle(fromY(e), 'quiet'); });
const tu = () => { this.thrDrag = false; $('#gl').focus({ preventScroll: true }); };
tr.addEventListener('pointerup', tu); tr.addEventListener('pointercancel', tu);
this.bindInput();
},
renderHangar() {
const sh = this.ship;
$('#hg-list').innerHTML = SHIPS.map((s, i) => `<button role="radio" aria-checked="${i === this.shipIdx}" class="hg-ship ${i === this.shipIdx ? 'on' : ''}" data-i="${i}"><b>${s.name}<span class="vb" data-vbs="${s.id}">${UI.shipBadge(s)}</span></b><small>${s.type}${s.info && s.info.franquicia ? ' · ' + s.info.franquicia : ''}</small></button>`).join('');
document.querySelectorAll('.hg-ship').forEach(b => b.addEventListener('click', () => { this.setShip(+b.dataset.i); this.renderHangar(); }));
$('#hg-type').textContent = sh.type; $('#hg-name').textContent = sh.name; $('#hg-desc').textContent = sh.desc;
const bar = (v, max, log) => { const k = log ? Math.log10(v) / Math.log10(max) : v / max; return `<span class="bar"><i style="transform:scaleX(${clamp(k, 0.04, 1)})"></i></span>`; };
const maneu = sh.angAcc >= 80 ? 'alta' : sh.angAcc >= 30 ? 'media' : 'baja';
$('#hg-stats').innerHTML = `
<dt>Velocidad máxima</dt><dd>${fmt(sh.maxSpeed, 0)} km/s${sh.maxSpeed > 3000 ? ' (' + fmt(sh.maxSpeed / C_KMS * 100, 1) + ' % de la luz)' : ''}${bar(sh.maxSpeed, 60000, true)}</dd>
<dt>Aceleración</dt><dd>${fmt(sh.accel, 0)} m/s² (${fmt(sh.accel / 9.81, 1)} g)${bar(sh.accel, 2500, true)}</dd>
<dt>Maniobrabilidad</dt><dd>${fmt(sh.angAcc, 0)} °/s² (${maneu})${bar(sh.angAcc, 110)}</dd>
<dt>Masa</dt><dd>${sh.massEst ? '≈ ' : ''}${fmt(sh.mass, 0)} t${sh.massEst ? ' <small class="est">estimación para la simulación</small>' : ''}</dd>
<dt>Longitud</dt><dd>${sh.lenEst ? '≈ ' : ''}${fmt(sh.len, 0)} m${sh.lenEst ? ' <small class="est">estimación</small>' : ''}</dd>
<dt>Sistema de propulsión</dt><dd>${sh.prop}</dd>
<dt>Viaje Tierra–Marte</dt><dd>≈ ${fmtDur(this.tripDays(sh, 0.52 * AU_KM))} de tiempo simulado</dd>`;
const I = sh.info, row = (k, v) => v ? `<dt>${k}</dt><dd>${esc(v)}</dd>` : '';
$('#hg-extra').innerHTML = (I ? `
<h4>Origen</h4><dl class="hg-dl">${row('Universo / serie', I.franquicia)}${row('Afiliación', I.afiliacion)}</dl>
<h4>Especificaciones</h4><dl class="hg-dl">${row('Clase / tipo', I.clase)}${row('Función', I.funcion)}${row('Tripulación', I.tripulacion)}${row('Velocidad', I.velocidad)}${row('Capacidad operativa', I.capacidad)}${row('Armamento', I.armamento)}</dl>
${I.destacadas ? `<h4>Características destacadas</h4><ul class="hg-feat">${I.destacadas.map(x => `<li>${esc(x)}</li>`).join('')}</ul>` : ''}` : '')
+ `<h4>Estado del modelo</h4><p class="hg-model" data-vbs-full="${sh.id}">${UI.shipBadge(sh, true)}</p>`;
const br = $('#hg-brand');
br.hidden = !sh.brand;
if (sh.brand) { const im = br.querySelector('img'); if (im.getAttribute('src') !== sh.brand.src) im.src = sh.brand.src; im.alt = sh.brand.alt; im.width = sh.brand.w; im.height = sh.brand.h; br.dataset.brand = sh.id; }
const cr = MODEL_CREDITS[sh.id]; $('#hg-credit').innerHTML = cr && World.hd[sh.id] ? 'Modelo 3D: ' + UI.creditHTML(cr) : '';
$('#hg-start').value = this.startId;
document.querySelectorAll('[data-as]').forEach(b => { const on = (b.dataset.as === '1') === this.assisted; b.classList.toggle('on', on); b.setAttribute('aria-checked', on); });
},
tripDays(sh, km) {
const a = sh.accel / 1000, v = sh.maxSpeed, dAcc = v * v / a;
const t = km < dAcc ? 2 * Math.sqrt(km / a) : 2 * v / a + (km - dAcc) / v;
return t / 86400;
},
toggleTargets(v) {
const p = $('#fl-targets'), open = v == null ? !p.classList.contains('open') : v;
p.classList.toggle('open', open);
if (open) { $('#ft-q').value = ''; this.renderTargets(); $('#ft-q').focus(); } else $('#gl').focus({ preventScroll: true });
},
renderTargets() {
const q = norm($('#ft-q').value.trim());
const groups = [['Planetas y el Sol', r => r.isSun || r.def.type === 'planet'], ['Lunas', r => r.def.type === 'moon'], ['Planetas enanos y transneptunianos', r => r.def.type === 'dwarf' || r.def.type === 'tno'], ['Asteroides', r => r.def.type === 'asteroid'], ['Cometas', r => r.def.type === 'comet'], ['Estaciones y naves', r => r.isCraft && !r.hidden]];
const match = r => !q || [r.def.name, r.def.short || '', ...(r.def.aka || [])].some(n => norm(n).includes(q));
$('#ft-list').innerHTML = groups.map(([t, f]) => {
const items = World.rb.filter(r => f(r) && match(r));
if (!items.length) return '';
return `<h5>${t}</h5>` + items.map(r => `<button data-t="${r.id}" class="${r === this.target ? 'on' : ''}"><i class="dot" style="--c:${r.def.color}"></i><span>${esc(r.def.short || r.def.name)}${UI.badge(r)}</span><small>${fmtKm(V.dist(this.pos, r.km))}</small></button>`).join('');
}).join('') || '<p class="note">Sin resultados.</p>';
document.querySelectorAll('#ft-list [data-t]').forEach(b => b.addEventListener('click', () => { this.setTarget(World.byId[b.dataset.t]); this.toggleTargets(false); }));
},
refreshBar() {
if (!this.dom) return;
const btn = (id, label, on, tip, key) => `<button class="fb ${on ? 'on' : ''}" data-fb="${id}" aria-pressed="${!!on}" data-tip="${tip} (${key})">${label}</button>`;
$('#fl-bar').innerHTML =
`<span class="fb-grp">` +
`<button class="fb cam" data-fb="cam" data-tip="Cambiar cámara (C)">${CAMS[this.camIdx][1]}</button>` +
btn('mode', this.assisted ? 'Asistido' : 'Simulación', !this.assisted, 'Modo de vuelo', 'V') +
btn('stab', 'Estabilizar', this.stab, 'Estabilización automática', 'Z') +
`</span><span class="fb-grp">` +
`<button class="fb tgt" data-fb="tgt" data-tip="Elegir objetivo (T)">${this.target ? esc(this.target.def.short || this.target.def.name) : 'Sin objetivo'}</button>` +
btn('aim', 'Orientar', this.aim, 'Orientar hacia el objetivo', 'G') + btn('match', 'Igualar', this.match, 'Igualar velocidad con el objetivo', 'B') +
btn('intercept', 'Interceptar', this.intercept, 'Piloto automático hasta el objetivo', 'I') + btn('hold', 'Mantener', !!this.hold, 'Mantener la distancia actual', 'H') +
btn('orbit', 'Órbita', !!this.orbit, 'Entrar o salir de órbita', 'O') + btn('autoBrake', 'Freno aprox.', this.autoBrake, 'Freno automático de aproximación', 'N') +
`</span><span class="fb-grp">` +
`<span class="fb-warpctl" role="group" aria-label="Aceleración temporal">` +
`<button class="fb fb-step" data-fb="w-" data-tip="Menos aceleración temporal ( , )" aria-label="Menos aceleración temporal">−</button>` +
`<span class="fb-warp" data-tip="Aceleración temporal" aria-live="polite">×${fmt(WARPS[this.warpIdx], 0)}</span>` +
`<button class="fb fb-step" data-fb="w+" data-tip="Más aceleración temporal ( . )" aria-label="Más aceleración temporal">+</button>` +
`</span>` +
btn('snd', FlightAudio.on ? 'Sonido' : 'Silencio', FlightAudio.on, 'Sonido de cabina', 'M') +
`<button class="fb ${Music.prefs.on && (Music.playing || Music.waiting) ? 'on' : ''}" data-fb="music" aria-pressed="${!!(Music.prefs.on && (Music.playing || Music.waiting))}" data-tip="Música de fondo: ${MUSIC_INFO.title}">Música</button>` +
`<button class="fb" data-fb="pause" data-tip="Pausa y ayuda (Esc)">Pausa</button>` +
`<span class="fb-touch"><button class="fb hold" data-hold="Space">Frenar</button></span>` +
`</span>`;
document.querySelectorAll('[data-fb]').forEach(b => b.addEventListener('click', () => {
const a = b.dataset.fb;
if (a === 'cam') this.cycleCam(1); else if (a === 'mode') this.toggleAssisted(); else if (a === 'tgt') this.toggleTargets();
else if (a === 'w-') this.setWarp(this.warpIdx - 1); else if (a === 'w+') this.setWarp(this.warpIdx + 1);
else if (a === 'snd') this.toggleSound(); else if (a === 'music') Music.toggle(); else if (a === 'pause') this.pause(true);
else this.toggle(a);
$('#gl').focus({ preventScroll: true });
}));
document.querySelectorAll('[data-hold]').forEach(b => {
const k = b.dataset.hold === 'Space' ? Keys.code('flight', 'brake') : b.dataset.hold, on = e => { e.preventDefault(); this.keys.add(k); }, off = () => this.keys.delete(k);
b.addEventListener('pointerdown', on); b.addEventListener('pointerup', off); b.addEventListener('pointerleave', off); b.addEventListener('pointercancel', off);
});
},
fmtSpeed(v) {
if (v < 1) return fmt(v * 1000, v < 0.1 ? 1 : 0) + ' m/s';
return fmt(v, v < 100 ? 2 : 0) + ' km/s' + (v > C_KMS * 0.001 ? ' (' + fmt(v / C_KMS * 100, 2) + ' % c)' : '');
},
hud() {
if (this.state === 'hangar' || !this.dom) return;
const W = innerWidth, H = innerHeight, ref = this.ref, T = this.target;
const vrel = V.sub(this.vel, ref.vel), sp = V.len(vrel);
const put = (el, x, y, on) => { if (!on || !isFinite(x) || !isFinite(y)) { el.classList.remove('on'); return; } el.style.transform = `translate(${x.toFixed(1)}px, ${y.toFixed(1)}px)`; el.classList.add('on'); };
const projDir = (d, o) => World.project(V.add(World.cam, d), o);
const o = {};
projDir(this.F, o); put($('#mk-bore'), o.x, o.y, o.on && this.camMode !== 'free');
if (sp > 1e-5) { const vd = V.scale(vrel, 1 / sp); projDir(vd, o); put($('#mk-pro'), o.x, o.y, o.on); projDir(V.scale(vd, -1), o); put($('#mk-retro'), o.x, o.y, o.on); }
else { $('#mk-pro').classList.remove('on'); $('#mk-retro').classList.remove('on'); }
const mt = $('#mk-tgt'), ar = $('#mk-arrow');
if (T) {
const Ts = T.posS; World.project(Ts, o);
const d = V.dist(this.pos, T.km), rv = V.sub(this.vel, T.vel), rh = V.scale(V.sub(T.km, this.pos), 1 / Math.max(d, 1e-9)), clos = V.dot(rv, rh);
const inView = o.on && o.x > 30 && o.x < W - 30 && o.y > 30 && o.y < H - 30;
if (inView) {
const s = Math.max(28, Math.min(T.proj.rpx * 2 + 18, 260));
mt.style.width = mt.style.height = s + 'px';
put(mt, o.x - s / 2, o.y - s / 2, true); ar.classList.remove('on');
if (this.frameN % 3 === 0) { mt.querySelector('b').textContent = T.def.short || T.def.name; mt.querySelector('span').textContent = fmtKm(d) + '  ' + (clos >= 0 ? '−' : '+') + this.fmtSpeed(Math.abs(clos)); }
} else {
mt.classList.remove('on');
const v = World.view, dir = V.sub(Ts, World.cam);
const vx = v[0] * dir[0] + v[4] * dir[1] + v[8] * dir[2], vy = v[1] * dir[0] + v[5] * dir[1] + v[9] * dir[2];
const ang = Math.atan2(-vy, vx), cx = W / 2, cy = H / 2, rx = W / 2 - 60, ry = H / 2 - 70;
const k = Math.min(rx / Math.abs(Math.cos(ang) || 1e-9), ry / Math.abs(Math.sin(ang) || 1e-9));
ar.style.transform = `translate(${(cx + Math.cos(ang) * k).toFixed(1)}px, ${(cy + Math.sin(ang) * k).toFixed(1)}px) rotate(${ang}rad)`;
ar.classList.add('on');
}
} else { mt.classList.remove('on'); ar.classList.remove('on'); }
this.updateThrottleUI();
this.frameN++;
if (this.frameN % 4 && !this.forceHud) return this.hudWarn();
this.forceHud = false;
const dl = (rows) => rows.map(([k, v, s]) => `<div><dt>${k}</dt><dd>${v}${s ? `<small>${s}</small>` : ''}</dd></div>`).join('');
const north = [0, 1, 0], F = this.F;
const z0 = x => (Math.abs(x) < 0.5 ? 0 : x);
const hdg = (vv) => { const e = Astro.sceneToEcl(vv); return wrap360(Math.round(Math.atan2(e[1], e[0]) * RAD)) % 360; };
const elev = (vv) => z0(Math.asin(clamp(V.dot(V.norm(vv), north), -1, 1)) * RAD);
const rollRef = V.norm(V.sub(north, V.scale(F, V.dot(north, F))));
const roll = z0(Math.atan2(V.dot(V.cross(rollRef, this.U), F), V.dot(rollRef, this.U)) * RAD);
const drift = sp > 1e-6 ? Math.acos(clamp(V.dot(V.norm(vrel), F), -1, 1)) * RAD : 0;
const E = World.byId.tierra;
$('#h-nav').innerHTML = dl([
['Velocidad', this.fmtSpeed(sp), 'respecto a ' + ref.def.name],
['Dirección de movimiento', sp > 1e-6 ? 'Rumbo ' + fmt(hdg(vrel), 0) + '°, elevación ' + fmt(elev(vrel), 0) + '°' : 'Detenida', sp > 1e-6 ? fmt(drift, 0) + '° respecto a la proa' : ''],
['Orientación', 'Rumbo ' + fmt(hdg(F), 0) + '°, cabeceo ' + fmt(elev(F), 0) + '°', 'Alabeo ' + fmt(roll, 0) + '° (marco eclíptico)'],
['Distancia al Sol', fmtKm(V.len(this.pos)), fmtAU(V.len(this.pos) / AU_KM)],
['Distancia a la Tierra', fmtKm(V.dist(this.pos, E.km)), V.dist(this.pos, E.km) > 1e6 ? fmtAU(V.dist(this.pos, E.km) / AU_KM) : ''],
]);
if (T) {
const d = V.dist(this.pos, T.km), rv = V.sub(this.vel, T.vel), rh = V.scale(V.sub(T.km, this.pos), 1 / Math.max(d, 1e-9)), clos = V.dot(rv, rh);
const alt = d - T.R, eta = clos > 1e-6 ? alt / clos : null;
$('#h-tgt').innerHTML = dl([
['Objetivo', esc(T.def.name), TYPE_LABEL[T.def.type]],
['Distancia', fmtKm(d), (T.R > 1 ? 'Altitud ' + fmtKm(Math.max(alt, 0)) : '') + (d > 1e6 ? '  ' + fmtAU(d / AU_KM) : '')],
['Velocidad relativa', this.fmtSpeed(V.len(rv)), clos >= 0 ? 'Acercándose a ' + this.fmtSpeed(clos) : 'Alejándose a ' + this.fmtSpeed(-clos)],
['Tiempo estimado de llegada', eta != null ? fmtDur(eta / 86400) : 'Sin aproximación', eta != null && WARPS[this.warpIdx] > 1 ? '≈ ' + fmtDur(eta / 86400 / WARPS[this.warpIdx]) + ' reales a ×' + fmt(WARPS[this.warpIdx], 0) : ''],
['Distancia al Sol', fmtKm(V.len(T.km)), fmtAU(V.len(T.km) / AU_KM)],
]);
} else $('#h-tgt').innerHTML = '<p class="note">Pulsa T o haz clic sobre un planeta, luna o estación para fijarlo como destino.</p>';
$('#h-ship').textContent = this.ship.name;
$('#h-mode').textContent = this.assisted ? 'Modo asistido' : 'Modo simulación';
$('#h-cam').textContent = 'Cámara ' + CAMS[this.camIdx][1].toLowerCase();
$('#h-warp').textContent = 'Tiempo ×' + fmt(WARPS[this.warpIdx], 0) + (this.maxWarpIdx < WARPS.length - 1 ? ' (máx. ×' + fmt(WARPS[this.maxWarpIdx], 0) + ')' : '');
this.hudWarn();
},
updateThrottleUI() {
const auto = this.intercept || this.match || this.hold || this.orbit;
const act = auto ? this.thrAuto : this.thrAct, cmd = auto ? null : this.thrCmd;
const pos = v => (v + 0.25) / 1.25 * 100, zero = pos(0);
const fill = $('#thr-fill');
fill.style.bottom = (act >= 0 ? zero : pos(act)) + '%'; fill.style.height = Math.abs(act) / 1.25 * 100 + '%';
fill.classList.toggle('rev', act < 0);
const lv = $('#thr-lever'); lv.style.bottom = pos(cmd == null ? act : cmd) + '%'; lv.classList.toggle('auto', !!auto);
const pa = Math.round(act * 100), pc = cmd == null ? null : Math.round(cmd * 100);
if (pa !== this._pa) { $('#thr-act').textContent = (pa < 0 ? '−' + -pa : pa) + ' %'; this._pa = pa; }
const ct = auto ? 'Piloto automático' : (pc < 0 ? '−' + -pc : pc) + ' %';
if (ct !== this._pc) { $('#thr-cmd').textContent = ct; $('#thr-track').setAttribute('aria-valuenow', pc == null ? pa : pc); this._pc = ct; }
if (this.frameN % 4 === 0) {
$('#thr-acc').textContent = fmt(this.accNow, this.accNow < 10 ? 2 : 0) + ' m/s² (' + fmt(this.accNow / 9.81, 2) + ' g)';
$('#thr-name').textContent = auto ? '' : (DETENT_NAMES[String(this.thrCmd)] || '');
}
},
hudWarn() {
const now = performance.now();
this.msgs = this.msgs.filter(m => m.until > now);
const items = [];
if (this.alert) items.push(`<p class="w alert">${esc(this.alert.txt)}</p>`);
if (this.braking) items.push('<p class="w">Freno de aproximación activo</p>');
if (this.warpLimited > now) items.push('<p class="w">Aceleración temporal reducida por proximidad</p>');
if (this.orbit) items.push(`<p class="w ok">En órbita de ${esc(this.orbit.B.def.name)} a ${fmtKm(this.orbit.r - this.orbit.B.R)} de altitud</p>`);
if (this.intercept) items.push('<p class="w ok">Piloto automático: interceptando</p>');
this.msgs.slice(-3).forEach(m => items.push(`<p class="w ${m.kind === 'alert' ? 'alert' : m.kind === 'warn' ? '' : 'ok'}">${esc(m.t)}</p>`));
const html = items.join('');
if (html !== this.lastWarn) { $('#hud-warn').innerHTML = html; this.lastWarn = html; }
},
};
})();