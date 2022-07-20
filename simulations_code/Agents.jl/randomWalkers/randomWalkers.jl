using Agents
using GLMakie

#Model constants
const repulsion_force = 1
const diameter = 1.
const dt = .1
const β = 0.1
const extent = (12.,12.)
const n_walkers = 20
const n_steps = 1000
const D = 1.
const save_step = 100

#Create agent
mutable struct randomWalker <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    f::NTuple{2,Float64}
end

#Initialization
function initialize_model(n_walkers)
    space2d = ContinuousSpace(extent, 3*diameter)

    model = ABM(randomWalker, space2d)

    for _ in 1:n_walkers
        pos = (extent[1]*rand(),extent[2]*rand())
        vel = (0,0)
        f = (0,0)

        add_agent!(pos,model,vel,f)
    end

    return model
end

#Compute forces
function compute_forces!(walker, model)

    neighbor_ids = nearby_ids(walker,model,diameter)

    F = (0.,0.)
    for id in neighbor_ids

        #Compute distance between agents
        neighbor = model[id].pos
        d = sqrt(sum((walker.pos .-neighbor).^2))

        #Add to repulsion force if closeby
        if d < diameter

            v = neighbor .- model[id].pos
            modulus = sqrt(sum(v.^2))
            if modulus > 0
                F = F .- repulsion_force .* v ./modulus
            end

        end

        walker.f = F

    end
end

#Integrator
function euler!(walker)

    # Make an euclidean step
    walker.pos = walker.pos .+ walker.vel .* dt
    walker.vel = Tuple(walker.vel .+ walker.f .*dt .- β .* walker.vel .+ sqrt(2D) .* randn(2) .*sqrt(dt))

    #Check boundaries
    if walker.pos[1] < 0
        walker.pos = (0.1,walker.pos[2])
    elseif walker.pos[1] > extent[1]
        walker.pos = (extent[1]-.1,walker.pos[2])
    end

    if walker.pos[2] < 0
        walker.pos = (walker.pos[1],0.1)
    elseif walker.pos[2] > extent[2]
        walker.pos = (walker.pos[1],extent[2]-.1)
    end

end

#step
function step!(model)
    for id in Schedulers.fastest(model)
        compute_forces!(model[id],model)
    end
    for id in Schedulers.fastest(model)
        euler!(model[id])
    end
end

#Run model without saving
model = initialize_model(n_walkers)
agent, model = run!(model,dummystep,step!,n_steps,adata=[:pos],when=range(1,n_steps,step=save_step))

#Plots
for step in unique(agent.step)
    x = [i[1] for i in agent[agent.step .== step,:].pos]
    y = [i[2] for i in agent[agent.step .== step,:].pos]
    fig = Figure(figsize=[10,10])
    ax = Axis(fig[1,1])
    meshscatter!(ax,x,y,markersize=diameter/2)
    xlims!(0,extent[1])
    ylims!(0,extent[2])

    save(string("outcome/",step,".png"),fig)
end

