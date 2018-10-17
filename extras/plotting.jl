using PGFPlotsX
using LaTeXStrings
PGFPlotsX.latexengine!(PGFPlotsX.PDFLATEX)

function timeplot(t, x; ylabel)
    @pgf Axis(
        {
            xmajorgrids,
            ymajorgrids,
            enlargelimits=false,
            enlarge_y_limits,
            xlabel=L"$t$ [s]",
            ylabel=ylabel
        },
        Plot(Table([t, x]))
    )
end

function frac_to_latex(frac::Rational)
    absfrac = abs(frac)
    ret = LaTeXString("\\frac{$(absfrac.num)}{$(absfrac.den)}")
    if frac < 0
        ret = "-" * ret
    end
    ret
end

function piticks!(p, pifrac::Rational; axis::Symbol)
    pifracrange = UnitRange(round.(Int, extrema(θ) ./ (pifrac * π), (RoundDown, RoundUp))...)
    labels = map(pifracrange) do i
        if i == 0
            L"$0$"
        else
            frac = pifrac * i
            frac == 1 ? L"$\pi$" : latexstring(frac_to_latex(pifrac * i) * "\\pi")
        end
    end
    ticks = pifracrange * pifrac * π
    p[(x = "xtick", y = "ytick")[axis]] = ticks
    p[(x = "xticklabels", y = "yticklabels")[axis]] = labels
    p
end
