# SOLARIS

Simulador educativo en 3D del Sistema Solar. Es un sitio estático (`index.html`, `css/solaris.css`, `js/solaris.js` y `assets/`) que GitHub Pages publica directamente desde la rama `main`: cada merge a `main` se despliega solo.

## Versión del sitio

Cada cambio que afecte al sitio publicado (HTML, CSS, JS o recursos de `assets/`) debe actualizar la versión de SOLARIS en el mismo commit:

- La versión está en `js/solaris.js`, en `const SOLARIS_BUILD = '…';`. El sitio la muestra en Ajustes › General ("Versión de SOLARIS: …") y en la consola del navegador.
- Formato: `AAAA.MM.DD-HHMM-web`, con la fecha y hora **UTC** del momento del cambio. Se obtiene con:

  ```sh
  date -u +%Y.%m.%d-%H%M-web
  ```

- Los cambios que no tocan el sitio (por ejemplo, solo este archivo) no necesitan nueva versión.

## Parámetros de caché

`index.html` carga la hoja de estilos y el script con un parámetro `?v=` que son los primeros 10 caracteres del md5 del archivo. Tras modificar `css/solaris.css` o `js/solaris.js` (incluido el cambio de versión), hay que recalcularlo para que los navegadores descarguen el archivo nuevo:

```sh
C=$(md5sum css/solaris.css | cut -c1-10); J=$(md5sum js/solaris.js | cut -c1-10)
sed -i "s#css/solaris.css?v=[0-9a-f]*#css/solaris.css?v=$C#; s#js/solaris.js?v=[0-9a-f]*#js/solaris.js?v=$J#" index.html
```

Actualiza la versión **antes** de recalcular el `?v=` del script, porque cambiar la versión cambia el md5 de `js/solaris.js`.
