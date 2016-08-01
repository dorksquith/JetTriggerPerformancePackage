echo "\documentclass{beamer}
\usepackage[english]{babel}
\usepackage{subfigure}

\begin{document}

\title{Results Jet Trigger Performance Package}
%\author{Edgar Kellermann for the Trigger-Level Analysis team}
\author{Edgar Kellermann}
\date{\today}

\begin{frame}
  \titlepage
\end{frame}" > results.tex

for var in "$@"
do
echo "\begin{frame}
    \begin{figure}[h]
      \centering
      \includegraphics[scale=0.50]{"$var"}
      \caption{"$var"}
      \label{img:1}
    \end{figure}
\end{frame}" >> results.tex
done

echo "\end{document}" >> results.tex

pdflatex results.tex
pdflatex results.tex


