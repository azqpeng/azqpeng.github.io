\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}
\usepackage{biblatex}
\usepackage{graphicx}
\graphicspath{ {./pluckers.PNG/} }

\title{Mathematical music}
\author{Alexander Peng}
\date{December 2022}

\begin{document}

\maketitle

In this report, we'll engage in two separate studies on the mathematics of music. One attempts to formally compute how certain musical instruments work, while the other takes a broader conceptual swipe at music theory. I hope you'll enjoy them both.

\section{How do you pluck a string?}

What types of waves are produced by plucking a string on a musical instrument? To answer this question, we'll need to figure out what plucking a string looks like, and how to break it down into recognizable wave patterns. One will assume some prior knowledge of Fourier series and differential equations, as it's fairly difficult approaching this problem without those tools.
\begin{figure}[h]
\caption{Depiction of pluck}
\centering
\includegraphics[scale=0.3]{plucker}
\end{figure}
\\
First, let's define our problem. To do this calculation, we'll assume the instrument is a simple tense string and that it is played by being plucked - like a single guitar string. The string itself has some arbitrary tension $T$, linear mass density $\mu$, length $L$ with two fixed endpoints. Assume that the string is plucked at distance $\frac{L}{2}$ from the origin with  some height $h$. For simplicity, let's set $L=1$. We'll abstractly approach the initial plucked string as a triangular wave, with the vertices of the triangle being the fixed endpoints and the point of pluck. \\

\subsection{Initial setup}

Using the given information we define the boundary conditions and initial conditions of the function $y(x,t)$ representing the displacement of the string as a function of time and location on the string. Because the string is fixed at both ends, displacement must always be zero at the ends.\\
First, the domain of $x$ is the closed  interval on the real line $[0,1]$. When $t=0$, we can observe that the following is true:\\
\begin{equation}
    \left\{
        \begin{aligned}
            & y(x,0) = \frac{2hx}{L}
            \quad 0 \leq x \leq 0.5\\
            & y(x,0) = 2h(1 - \frac{x}{L}) \quad 0.5 < x \leq 1
        \end{aligned} 
    \right.
\end{equation}
Second, if we assume that the pluck is perfect, it should be fair to assume that the initial velocity is zero, or:\\
\begin{equation}
    \frac{dy(x,t)}{dt} |_{t=0} = 0    
\end{equation}
Finally, this is a string, giving us a transverse string wave, which means we can define it with the wave equation: 
\begin{equation}
    \frac{\partial^2y}{\partial t ^2} = \frac{1}{v^2}\frac{\partial^2y}{\partial x ^2} 
\end{equation}


\subsection{Initial fourier series}

 Note that the function $y(x,t)$ is both periodic due to it's boundary conditions, as $0 = y(0, t) = y(L, t)$ and continuous (though... it's not a particularly differentiable function). Thus, one can apply the fourier theorem to decompose it into a countably infinite sum of sinusoidal waves. Let's suppose we set $t=0$ to find the fourier series for initial state. Applying the formulation directly, we get:\\
\begin{equation*}
    \begin{aligned}
        y(x,0) = \sum_{n=-\infty}^{\infty} c_n e^{\frac{i\pi nx}{L}} \quad \forall n \in \mathbf{Z}\\
    \end{aligned}
\end{equation*}
One can note here that the wave obtained, $y$, is a triangle wave with maximum center height $d$, meaning one half of the wave has slope $\frac{2d}{L}$ and the other has slope $- \frac{2d}{L}$. This is trivially equivalent to the integral of a square wave of amplitude $\pm\frac{2d}{L}$ over the defined domain, cutting off in the middle, which we can define as such:\\
\begin{equation}
    \left\{
        \begin{aligned}
            & sq(x,0) = \frac{2h}{L} \quad 0 \leq x \leq 0.5\\
            & sq(x,0) = -\frac{2h}{L} \quad 0.5 < x \leq 1
        \end{aligned} 
    \right.
\end{equation}
If one expresses this as a Fourier series, then:\\
\begin{equation}
        \begin{aligned}
        c_k = \int_0^{0.5} 2he^{-i2\pi nt} + \int_{0.5}^1 -2he^{-i2\pi nt}
        = -\frac{4h}{i2\pi n}((-1)^n - 1)\\
        \end{aligned} 
\end{equation}
Or equivalently:\\
\begin{equation}
    \left\{
        \begin{aligned}
            & \frac{4h}{i\pi n} \quad \text{if n odd}\\
            & 0 \quad \text{otherwise}
        \end{aligned}
    \right.
\end{equation}
Putting it all together by taking the integral of this Fourier expression and doubling the period, as waves have both positive and negative components, one obtains:\\
\begin{equation}
    y(x,0) = \int_{-\infty}^{x} \sum_{\text{$n$ odd}} \frac{4h}{i\pi n} e^{i\pi nx} = \sum_{\text{$n$ odd}} -\frac{4h}{\pi^2 n^2} e^{i\pi nx} = \sum_{n=0}^{\infty} \frac{8h}{\pi^2 n^2}sin(\frac{\pi n}{2})sin(\pi nx)
\end{equation}
Note that the last step is done by symmetry and using Euler's formula. Now, we have an expression, so we just need to put it together by adding in time (yay, waves)!

\subsection{Put it together}

Notice that there is no currently time dependence on the wave. While there are  tedious ways to calculate it, one can use physics to think of a simpler solution by solving the wave equations. Since no part of $y(x,t)$ that we've derived so far contains time dependencies, and plucks form vibrations (which are thus time/space independent), let's arbitrarily add a function $f(t)$ to account for the temporal component and simply solve the wave equation which we know must hold since we're dealing with a transverse string wave.
\begin{equation}
    \begin{aligned}
    & \frac{\partial^2y}{\partial t ^2} = \frac{1}{v^2}\frac{\partial^2y}{\partial x ^2}\\
    & \rightarrow \frac{\partial^2 \sum_{n=0}^{\infty} \frac{8h}{\pi^2 n^2}sin(\frac{\pi n}{2})sin(\pi nx)f(t)}{\partial x^2} = \frac{1}{v^2}\frac{\partial^2\sum_{n=0}^{\infty} \frac{8h}{\pi^2 n^2}sin(\frac{\pi n}{2})sin(\pi nx)f(t)}{\partial t^2}\\
    & \leftrightarrow \sum_{n=0}^{\infty} -8h sin(\frac{\pi n}{2})sin(\pi nx)f(t) = \sum_{n=0}^{\infty} \frac{8h}{v^2\pi^2 n^2}sin(\frac{\pi n}{2})sin(\pi nx)f(t)\\
    & \leftrightarrow f(t) = \frac{f''(t)}{v^2\pi^2 n^2}\\
    & \rightarrow f(t) = cos(v\pi n t) = cos(\omega t)\\
    \end{aligned}
\end{equation}
This just so happens to match the physical description of a standing wave, which reinforces that our strategy is working. From the equation derived previously, the harmonics for a transversal string wave with fixed ends occur at $k = \frac{2\pi}{\lambda} = \frac{2\pi}{2L/n} = \frac{n\pi}{L}$ and thus $\omega = kv = \frac{n\pi}{L}\sqrt{\frac{T}{\mu}}$. Consequently, the elements of our series should look like $c_n sin(kx)cos(\omega t)$, and judging from our calculation earlier, they do! Great!\\
\\
Now, for the big finale, we take our initial Fourier series, add in the time dependence we calculated, and substitute in $k$ and $\omega$ to understand it as a standing wave. The result is:\\
\begin{equation*}
    \begin{aligned}
        y(x,t) = \sum_{n=0}^{\infty} \frac{8h}{\pi^2 n^2}sin(\frac{\pi n}{2})sin(\pi nx)cos(v\pi n t) = \boxed{\sum_{n=0}^{\infty} \frac{8h}{\pi^2 n^2}sin(\frac{\pi n}{2})sin(kx)cos(\omega t)}
    \end{aligned}
\end{equation*}

\section{What are scales?}

NOTE: This section is mostly just documenting relationships between music, physics, and math. It is almost completely conceptual and not connected to the computational work above - but there were so many interesting music theory facts I learned along the way so I hoped to share some of them with you all.\\
\\
It's easy to assume that musical notes, the 7-step major scale, the 12-step chromatic scale, etc. are purely social constructs - and for the most part, they are. However, there are also mathematical and physical justifications for why these patterns work.\\
\\
Musicians have long noticed that human ears aren't respondent to frequency but respondent to pitch - which operates logarithmically. In other words, doubling frequency will sound like the same note, but tripling a note's frequency and doubling a note's frequency won't sound like the same note. However, quadrupling a note's frequency and doubling a note's frequency will. While at first one would assume that the reasoning for this is purely biological, it seems recent evidence actually suggests that octaves (doubling frequency) might be derived not in our genetics but through culture. I'm not qualified to speak on this, but I linked the article below (Renken).\\
\\
Regardless of why pitch and frequency have an exponential relationship to many human ears, this relationship means that defining music theory is quite hard. We'll examine this with a case study.\\
\\
Let's take the musical note $C_3$, which is a standing wave with fundamental frequency $220$Hz. What happens if you take the second harmonic of said frequency? At $440$Hz you get the note $A_4$, or an octave - a full scale - up from $A_3$. If you double the frequency of every note in a piece, it sounds like a higher version of the same thing. What happens if you take the third harmonic? At $660$Hz, you don't have $A_5$ - instead, you have approximately an $E_5$. When you play the second and third harmonics of a note together ($A_5$ and $E_5$), you get a very pleasant sound - in music theory, this is a "major fifth". Similarly, when you play the fourth and fifth harmonics of a note together, you get what is considered the most "consonant" interval - a "major third." Finally, when you play the eighth and ninth harmonic together, you obtain a musical "step" - 7 of which form a scale. The sounds formed by combinations of harmonic frequencies are great. Similarly, one of the most "dissonant" sounds in music theory - the tritone - takes a note and plays it with a note with frequency $\sqrt{2}$ times greater. That isn't close to any harmonic, and also happens to sound terrible. While there's rationale for the appeal of harmonics - notably pattern matching and preference for integers - there aren't many definitive cognitive studies on it.\\
\\
An interesting angle to take is to expand from music theory to the real world, say a piano. In theory, major scales require $7$ steps to go up an entire octave. However, if each step is the difference between an eighth and ninth harmonic, then the octave you create with seven steps is actually $(9/8)^7 = ~2.027286...$ times the frequency of the initial note, as opposed to double it. This is a small difference, but it is certainly nontrivial. If one wants to go up a scale by another common definition - 12 semitones or "half steps" - which are composed of fifteenth and sixteenth harmonics, then the ratio between first note and last note is $(16/15)^12 = ~2.169425...$ which isn't even close to flat doubling!\\
\\
The ultimate cause of this is that once again, notes are defined by pitch and pitch is exponential. And the $n$th roots of certain rational numbers, for example $2$, often are not rational, making it hard to be defined in terms of harmonics. To explore this further requires measure-theoretic understandings (punny!) of the real line, which lie outside the scope of this class - but a good introduction to those ideas is offered in the referencfes (Sanderson).\\
\\

    \section*{References}
    I like to think to myself why I picked a topic and how that informs the work I do on it - even if it's just a single homework assignment. I'm a math and computer science major who wants to apply the types of math I like (in this case, real/harmonic analysis and differential equations) in real world settings. So I wanted to investigate fourier series, and that was my main motivation. I used a decent amount of outside mathematical knowledge to do this, including things like dirac delta functions and of course fourier series. A secondary motivation is that, being interested in fourier analysis, I'm also interested in mathematical analysis writ large - and so wanted to do some things related to it, which I was hoping to do towards the end with the scale stuff but didn't really flesh out beyond conceptual ideas since I felt like I was drifting a bit out of the scientific world. Below are some sources I used:\\
    \\
    Johnson, Don. Fundamentals of electrical engineering I. OpenStax CNX, 2015.\\
    \\
    Sanderson, Grant. “Music and Measure Theory.” 3Blue1Brown, 3 Oct. 2015,\\ https://www.youtube.com/watch?v=cyW5z-M2yzw. \\
    \\
    Renken, Elena. “Perceptions of Musical Octaves Are Learned, Not Wired in the Brain.” Quanta Magazine, 5 June 2020, https://www.quantamagazine.org/perceptions-of-musical-octaves-are-learned-not-wired-in-the-brain-20191030/. \\
    \\
    “Intervals.” AudioLabs, International Audio Laboratories Erlangen, https://www.audiolabs-erlangen.de/resources/MIR/FMP/C5/C5S1_Intervals.html.\\
    \\
    minutephysics. “Why It's Impossible to Tune a Piano.” YouTube, YouTube, 17 Sept. 2015, https://www.youtube.com/watch?v=1Hqm0dYKUx4. 



\end{document}
