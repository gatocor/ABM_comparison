using AgentBasedModels
using GLMakie

#Model constants
const repulsion_force = 1
const diameter = 1
const dt = .1
const β = 0.1
const extent = (102.,102.,102.)
const n_walkers = 1
const tMax = 5000
const dtSave = 100
const div_rate = 0.001*dt

model = @agent(3,
    [div_rate,r,repulsion_force,β]::Global,
    [vx,vy,vz]::Local,
    [fx,fy,fz]::LocalInteraction,

    #Compute the pairwise forces
    UpdateInteraction=begin
        d = sqrt((x.i-x.j)^2+(y.i-y.j)^2+(z.i-z.j)^2)

        if 0 < d < 2*r

            dx = x.j .- x.i
            dy = y.j .- y.i
            dz = z.j .- z.i
            fx.i = fx.i .- repulsion_force .* dx ./d
            fy.i = fy.i .- repulsion_force .* dy ./d
            fz.i = fz.i .- repulsion_force .* dz ./d

        end
    end,

    #Differential equations describing the model
    UpdateVariable=begin
        d(x) = fx*dt - β*vx*dt
        d(y) = fy*dt - β*vy*dt
        d(z) = fz*dt - β*vz*dt
        d(x) = vx*dt
        d(y) = vy*dt
        d(z) = vz*dt
    end,

    #Make the division update rule
    UpdateLocal=begin
        if rand() < div_rate
            px = Normal(0,1)
            py = Normal(0,1)
            pz = Normal(0,1)
            mod = sqrt(px^2+py^2+pz^2)
    
            #New agent
            d = r^(1/3.)
            addAgent(x=x+d*px/mod,y=y+d*py/mod,z=z+d*pz/mod,vx=0,vy=0,vz=0)
            addAgent(x=x-d*px/mod,y=y-d*py/mod,z=z-d*pz/mod,vx=0,vy=0,vz=0)
            #Update original agent with other position
            removeAgent()
        end
    end
)
model = compile(model,platform="cpu",neighbors="full",integrator="Euler",save="RAM")

com = Community(model,N=1)
com.div_rate = div_rate
com.r = diameter/2
com.repulsion_force = repulsion_force
com.β = β
com.t = 0.

steps = Int64(round(tMax/dtSave))
for step in range(1,steps,step=1)

    @time comt = model.evolve(com, dt = dt, dtSave = dtSave, tMax = com.t+dtSave, nMax = 4000)
    global com = comt[end]
    println("Global Time: " , com.t , "/" , tMax)
    println("N agents: " , com.N , "\n")

end