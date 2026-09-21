# -*- coding: utf-8 -*-
r"""La cadena EXACTA del resumen para pegar en el formulario de arXiv, y su cuenta.

POR QUE EXISTE
--------------
arXiv corta el resumen en 1920 caracteres y este resumen tiene DIEZ de margen, asi que
la decision de como teclear una sola letra griega decide si el envio pasa o lo rechaza:

    griegas en Unicode (eta -> n griega, rho, gamma, Lambda)   1910  <- pasa
    griegas deletreadas ("eta", "rho", "gamma", "Lambda")      1940  <- RECHAZADO

docs/ARXIV_SUBMISSION.md mandaba a `scratchpad/cierre/abs_paste.py` para esto, que
vivia en el scratchpad de una sesion y ya no existe.  Este fichero lo repone dentro del
deposito, que es donde debia estar desde el principio: una instruccion de envio que
apunta a un script borrado no es una instruccion.

Y hay una trampa historica que este fichero NO repite.  `mide_resumen.py` borraba `^` y
`_` antes de contar, de modo que leia OCHO caracteres de menos: con ese contador el
resumen marcaba 1926 contra el limite de 1920 y se daba por bueno.  Aqui no se borra
nada antes de contar: se cuenta la cadena que se va a pegar, tal cual.

USO
---
    python src/abs_paste.py            # imprime las dos versiones y sus cuentas
    python src/abs_paste.py --write    # ademas deja build/abstract_arxiv.txt
"""
from __future__ import print_function

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
MAIN = os.path.join(REPO, "paper", "main.tex")
LIMITE = 1920

# --- griegas y simbolos, en Unicode: la forma que cabe
GRIEGAS = [
    (r"\varepsilon", u"\u03b5"), (r"\epsilon", u"\u03b5"),
    (r"\varphi", u"\u03c6"), (r"\phi", u"\u03c6"),
    (r"\eta", u"\u03b7"), (r"\rho", u"\u03c1"),
    (r"\gamma", u"\u03b3"), (r"\Lambda", u"\u039b"),
    (r"\lambda", u"\u03bb"), (r"\omega", u"\u03c9"),
    (r"\Omega", u"\u03a9"), (r"\chi", u"\u03c7"),
    (r"\sigma", u"\u03c3"), (r"\tau", u"\u03c4"),
    (r"\mu", u"\u03bc"), (r"\nu", u"\u03bd"),
    (r"\delta", u"\u03b4"), (r"\Delta", u"\u0394"),
    (r"\alpha", u"\u03b1"), (r"\beta", u"\u03b2"),
    (r"\theta", u"\u03b8"), (r"\pi", u"\u03c0"),
    (r"\psi", u"\u03c8"), (r"\Psi", u"\u03a8"),
]

# --- deletreadas: la forma que NO cabe, calculada para poder ensenar la diferencia
DELETREADAS = [(c, n.lstrip("\\")) for c, n in
               [(a, a) for a, _ in GRIEGAS]]

OTROS = [
    (r"\times", u"\u00d7"), (r"\approx", u"\u2248"),
    (r"\le", u"\u2264"), (r"\leq", u"\u2264"),
    (r"\ge", u"\u2265"), (r"\geq", u"\u2265"),
    (r"\ll", u"<<"), (r"\gg", u">>"),
    (r"\pm", u"\u00b1"), (r"\cdot", u"\u00b7"),
    (r"\dim", u"dim"), (r"\ldots", u"..."), (r"\dots", u"..."),
    (r"\,", u""), (r"\;", u" "), (r"\!", u""), (r"\ ", u" "),
    (r"\%", u"%"), (r"\&", u"&"), (r"\_", u"_"),
    (r"---", u"--"), (r"--", u"--"),
    (r"``", u'"'), (r"''", u'"'),
]


def extrae():
    t = io.open(MAIN, encoding="utf-8").read()
    m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", t, re.S)
    if not m:
        raise SystemExit("  no encuentro el entorno abstract en paper/main.tex")
    return m.group(1)


def a_texto(src, griegas_unicode=True):
    t = src

    # comentarios de LaTeX
    t = re.sub(r"(?m)^\s*%.*$", "", t)
    t = re.sub(r"(?<!\\)%.*", "", t)

    # ordenes con argumento cuyo CONTENIDO se conserva
    for c in ("emph", "textit", "textbf", "text", "mathrm", "mathcal", "mathbf",
              "operatorname", "mbox", "textrm"):
        for _ in range(4):
            t = re.sub(r"\\" + c + r"\s*\{([^{}]*)\}", r"\1", t)

    # exponentes y subindices: se CONSERVAN, porque cuentan
    t = re.sub(r"\^\{([^{}]*)\}", r"^\1", t)
    t = re.sub(r"_\{([^{}]*)\}", r"_\1", t)

    # fracciones simples
    t = re.sub(r"\\t?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}", r"\1/\2", t)

    if griegas_unicode:
        for c, u in GRIEGAS:
            t = t.replace(c + " ", u).replace(c, u)
    else:
        for c, _ in GRIEGAS:
            t = t.replace(c + " ", c.lstrip("\\") + " ").replace(
                c, c.lstrip("\\"))

    for c, u in OTROS:
        t = t.replace(c, u)

    # lo que quede de matematicas y agrupacion
    t = t.replace("$", "")
    t = re.sub(r"\\[a-zA-Z]+\*?", " ", t)
    t = t.replace("{", "").replace("}", "")
    t = t.replace("~", " ")

    # espacios
    t = re.sub(r"[ \t]+", " ", t)
    t = re.sub(r"\n\s*\n+", "\n\n", t)
    return t.strip()


def main():
    src = extrae()
    uni = a_texto(src, True)
    spel = a_texto(src, False)

    print("=" * 78)
    print("  RESUMEN PARA arXiv -- limite %d caracteres" % LIMITE)
    print("=" * 78)
    for nombre, s in (("griegas en Unicode", uni), ("griegas deletreadas", spel)):
        n = len(s)
        veredicto = "PASA, margen %d" % (LIMITE - n) if n <= LIMITE \
            else "RECHAZADO, sobran %d" % (n - LIMITE)
        print("  %-22s %5d caracteres   %s" % (nombre, n, veredicto))
    print()
    print("  palabras: %d   (APS pide menos de 500)" % len(uni.split()))
    print("  parrafos: %d   (APS pide UNO)" % len([p for p in uni.split("\n\n")
                                                   if p.strip()]))
    for mal, que in ((r"\begin{equation}", "ecuacion presentada"),
                     (r"\ref{", "referencia cruzada"),
                     (r"\cite{", "cita numerada")):
        if mal in src:
            print("  AVISO: el fuente contiene %s (%s): APS lo prohibe en el resumen"
                  % (mal, que))
    # El fichero SIEMPRE se escribe, y se escribe ANTES de imprimir.  La consola de
    # Windows es cp1252 y no sabe codificar una eta griega: la primera version de este
    # fichero imprimia primero y escribia despues, asi que la excepcion de consola se
    # llevaba por delante justo el fichero que hacia falta para el envio.
    d = os.path.join(REPO, "build")
    if not os.path.isdir(d):
        os.makedirs(d)
    p = os.path.join(d, "abstract_arxiv.txt")
    io.open(p, "w", encoding="utf-8", newline="\n").write(uni)
    print("  PEGAR DESDE: %s" % p)
    print("  (%d caracteres, UTF-8, sin BOM -- abrelo y copia todo)" % len(uni))
    print()
    print("-" * 78)
    try:
        print(uni)
    except UnicodeEncodeError:
        # la consola no puede; el fichero sí lo tiene bien
        print(uni.encode("ascii", "backslashreplace").decode("ascii"))
        print()
        print("  (la consola no codifica griegas: arriba van escapadas.")
        print("   El fichero de build/ las lleva correctas. Copia DEL FICHERO.)")
    print("-" * 78)

    return 0 if len(uni) <= LIMITE else 1


if __name__ == "__main__":
    sys.exit(main())
