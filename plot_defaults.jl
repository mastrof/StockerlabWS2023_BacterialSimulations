using Plots; Plots.plotlyjs()
Plots.default(
    thickness_scaling=1.5, grid=false, frame=:border,
    linewidth=1.5,
    guidefontsize=11, tickfontsize=11,
    bgcolor=:black, fgcolor=:white,
    fgcolorlegend=:transparent, bgcolorlegend=:transparent,
    palette=:Dark2,
    size=(750,500)
)