# line plot colors
red_plot = RGBf(233/255, 39/255, 56/255) 
grey_plot =  RGBf(206/255, 206/255, 206/255)
yellow_plot = RGBf(235/255, 158/255, 43/255) 

function plot_tracking_results((x_akmpc,u_akmpc,t_akmpc),(x_skmpc,u_skmpc,t_skmpc),  (x_lmpc,u_lmpc,t_lmpc), (x_ref, t_ref))
    
    @assert size(x_akmpc)[1] == size(x_skmpc)[1] == size(x_lmpc)[1] == size(x_ref)[1] "Dimension of input matrices corresponding to state dimensions need to be equal."
    linew = 3

    if size(x_akmpc)[1] == 2
        fig = Figure(fontsize = 30, size=(2200, 1200))
        ga = fig[1:2,1]= GridLayout() 
        
        # θ_1
        ax1b = Axis(ga[1,1], xlabel=L"t (s)",ylabel=L"$θ_1$ (rad)", width=1000,height=500,ytickformat = "{:.0f}")
        x1s = lines!(ax1b, t_skmpc[:], x_skmpc[1,:], color="gray", linewidth=linew)
        x1l = lines!(ax1b, t_lmpc[:], x_lmpc[1,:], color=yellow_plot, linewidth=linew/2)
        x1a = lines!(ax1b, t_akmpc[:], x_akmpc[1,:], color=red_plot, linewidth=linew)
        x1r = lines!(ax1b, t_ref[:], x_ref[1,:], color="black", linestyle=(:dash, :dense), linewidth=linew)
        
        # u_1 
        ax2b= Axis(ga[2,1], xlabel=L"t (s)",ylabel=L"$u_1$ (Nm)", width=1000,height=500,ytickformat = "{:.1f}")
        x1s = lines!(ax2b, t_skmpc[1:end-1], u_skmpc[1,:], color="gray", linewidth=linew)
        x1l = lines!(ax2b, t_lmpc[1:end-1], u_lmpc[1,:], color=yellow_plot, label=L"$$linearization",linewidth=linew)
        x1a = lines!(ax2b, t_akmpc[1:end-1], u_akmpc[1,:], color=red_plot, linewidth=linew)

        hidexdecorations!(ax1b, ticks = true, grid=false)
        rowgap!(ga, 0)
        colgap!(ga, 0)

        Legend(fig[1,2], [x1r, x1l, x1a, x1s], [L"$$reference", L"$$linearization MPC", L"$$adaptive KMPC", L"$$static KMPC"])
        # axislegend(ax1b, [x1r, x1l, x1a, x1s], [L"$$reference", L"$$linearization MPC", L"$$adaptive KMPC", L"$$static KMPC"], orientation = :horizontal)

    elseif size(x_akmpc)[1] == 4
        fig = Figure(fontsize = 30, size=(1100, 2100))
        ga = fig[1:2,1:2]= GridLayout() 
        ax1a = Axis(ga[1,1], xlabel=L"$t$ (s)",ylabel=L"$θ_2$ (rad)", width=1000,height=500,ytickformat = "{:.0f}")
        # θ_2
        x1a = lines!(ax1a, t_akmpc[:], x_akmpc[2,:], color=red_plot, linewidth=linew)
        x1l = lines!(ax1a, t_lmpc[:], x_lmpc[2,:], color=yellow_plot,linewidth=linew)
        x1s = lines!(ax1a, t_skmpc[:], x_skmpc[2,:], color="gray", linewidth=linew)
        x1r = lines!(ax1a, t_ref[:], x_ref[2,:], color="black", linestyle=(:dash, :dense), linewidth=linew)
        
        # θ_1
        ax1b = Axis(ga[2,1], xlabel=L"t (s)",ylabel=L"$θ_1$ (rad)", width=1000,height=500,ytickformat = "{:.0f}")
        x1a = lines!(ax1b, t_akmpc[:], x_akmpc[1,:], color=red_plot, linewidth=linew)
        x1l = lines!(ax1b, t_lmpc[:], x_lmpc[1,:], color=yellow_plot, linewidth=linew)
        x1s = lines!(ax1b, t_skmpc[:], x_skmpc[1,:], color="gray", linewidth=linew)
        x1r = lines!(ax1b, t_ref[:], x_ref[1,:], color="black", linestyle=(:dash, :dense), linewidth=linew)
        
        axislegend(ax1a, [x1r, x1l], [L"$$reference", L"$$linearization MPC"])
        axislegend(ax1b, [x1a, x1s], [L"$$adaptive KMPC", L"$$static KMPC"])
        hidexdecorations!(ax1a, ticks = true, grid=false)
    
        # u_2 
        ax2a= Axis(ga[1,2], xlabel=L"$t$ (s)",ylabel=L"$u_2$ (Nm)", width=1000,height=500,ytickformat = "{:.0f}")
        x1a = lines!(ax2a, t_akmpc[1:end-1], u_akmpc[2,:], color=red_plot, linewidth=linew)
        x1s = lines!(ax2a, t_skmpc[1:end-1], u_skmpc[2,:], color="gray",linewidth=linew)
        x1l = lines!(ax2a, t_lmpc[1:end-1], u_lmpc[2,:], color=yellow_plot,linewidth=linew)
        # u_1 
        ax2b= Axis(ga[2,2], xlabel=L"t (s)",ylabel=L"$u_1$ (Nm)", width=1000,height=500,ytickformat = "{:.0f}")
        x1a = lines!(ax2b, t_akmpc[1:end-1], u_akmpc[1,:], color=red_plot,linewidth=linew)
        x1s = lines!(ax2b, t_skmpc[1:end-1], u_skmpc[1,:], color="gray",linewidth=linew)
        x1l = lines!(ax2b, t_lmpc[1:end-1], u_lmpc[1,:], color=yellow_plot,linewidth=linew)
    
        hidexdecorations!(ax2a, ticks = true, grid=false)
        rowgap!(ga, 0)
    else    
        println("Dimension of input matrices corresponding to state dimensions can either be two or four.")
    end 

    resize_to_layout!(fig)
    return fig
end  

