module PlottingTools

export
    timeplot,
    piticks!

using PGFPlotsX
using LaTeXStrings
PGFPlotsX.latexengine!(PGFPlotsX.PDFLATEX)

function timeplot(t, x; ylabel, tmin=minimum(t), tmax=maximum(t))
    @pgf Axis(
        {
            xmajorgrids,
            ymajorgrids,
            xmin=tmin,
            xmax=tmax,
            enlarge_y_limits,
            xlabel=L"$t$ (s)",
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

function piticks!(p, pifrac::Rational, θ::AbstractVector; axis::Symbol)
    pifracrange = UnitRange(round.(Int, extrema(θ) ./ (pifrac * π), (RoundDown, RoundUp))...)
    labels = map(pifracrange) do i
        frac = pifrac * i
        if frac == 0
            L"$0$"
        elseif abs(frac) == 1
            if frac > 0
                L"$\pi$"
            else
                L"-$\pi$"
            end
        elseif frac.den == 1
            latexstring("$(frac.num) \\pi")
        else
            latexstring(frac_to_latex(pifrac * i) * "\\pi")
        end
    end
    ticks = pifracrange * pifrac * π
    p[(x = "xtick", y = "ytick")[axis]] = ticks
    p[(x = "xticklabels", y = "yticklabels")[axis]] = labels
    p
end

end # module
