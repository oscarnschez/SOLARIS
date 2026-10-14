# SOLARIS

Simulador educativo en 3D del Sistema Solar. Es un sitio estático (`index.html`, `css/`, `js/`, `assets/` y las páginas legales `terminos-y-condiciones/`, `politica-de-privacidad/` y `politica-de-cookies/`) que GitHub Pages publica directamente desde la rama `main`: cada merge a `main` se despliega solo.

## Versión del sitio

Cada cambio que afecte al sitio publicado (HTML, CSS, JS o recursos de `assets/`) debe actualizar la versión de SOLARIS en el mismo commit:

- La versión está en `js/solaris.js`, en `const SOLARIS_BUILD = '…';`. El sitio la muestra en Ajustes › General ("Versión de SOLARIS: …") y en la consola del navegador.
- Formato: `AAAA.MM.DD-HHMM-web`, con la fecha y hora **UTC** del momento del cambio. Se obtiene con:

  ```sh
  date -u +%Y.%m.%d-%H%M-web
  ```

- Los cambios que no tocan el sitio (por ejemplo, solo este archivo) no necesitan nueva versión.

## Parámetros de caché

`index.html` y las tres páginas legales cargan sus hojas de estilos y scripts con un parámetro `?v=` que son los primeros 10 caracteres del md5 del archivo. Tras modificar cualquiera de ellos (incluido el cambio de versión en `js/solaris.js`), hay que recalcularlos para que los navegadores descarguen el archivo nuevo:

```sh
for f in css/solaris.css css/fonts.css css/privacidad.css css/legal.css js/solaris.js js/privacidad.js; do
  h=$(md5sum "$f" | cut -c1-10)
  sed -i "s#$f?v=[0-9a-f]*#$f?v=$h#g" index.html terminos-y-condiciones/index.html politica-de-privacidad/index.html politica-de-cookies/index.html
done
```

Actualiza la versión **antes** de recalcular el `?v=` del script, porque cambiar la versión cambia el md5 de `js/solaris.js`.

## Privacidad y documentos legales

- SOLARIS no usa cookies, analítica ni servicios de terceros al cargar. Las tipografías se sirven desde `assets/fonts/` (`css/fonts.css`); no vuelvas a enlazar Google Fonts ni otros CDN.
- Todo contenido externo opcional (hoy, las fotografías de Wikimedia Commons) debe cargarse solo con permiso: `SolarisPrivacy.allows(...)` en `js/privacidad.js`, desactivado por defecto.
- Si añades una clave de `localStorage`, regístrala en `STORE` de `js/privacidad.js` y en la tabla de la Política de cookies.
- Si cambia el texto de un documento legal, sube su versión y fechas en la página y en `DOCS` de `js/privacidad.js`, y añade una fila a su historial de versiones. Los datos del responsable marcados como «Pendiente» no deben rellenarse con información supuesta.
