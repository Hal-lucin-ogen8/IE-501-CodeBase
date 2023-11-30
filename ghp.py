import mip

flights = ['xy', 'zy', 'wy', 'yw', 'xz', 'xw', 'yz', 'yx', 'zw', 'zx', 'wz', 'wx']
sectors = ['A', 'B', 'C', 'D', 'E', 'F']
airports = ['x', 'y', 'z', 'w']
lengths = {'x':1, 'y':1, 'z':1, 'w':1, 'A':1, 'B':2, 'C':1, 'D':2, 'E':3, 'F':4}

paths = {
    'xy': ['x', 'A', 'C', 'D', 'F', 'y'],
    'zy': ['z', 'B', 'D', 'F', 'y'],
    'wy': ['w', 'E', 'F', 'y'],
    'yw': ['y', 'F', 'E', 'w'],
    'xz': ['x', 'A', 'B', 'z'],
    'xw': ['x', 'A', 'C', 'D', 'E', 'w'],
    'yx': ['y', 'F', 'D', 'C', 'A', 'x'],
    'yz': ['y', 'F', 'D', 'B', 'w'],
    'zx': ['z', 'B', 'A', 'x'],
    'zw': ['z', 'B', 'E', 'w'],
    'wx': ['y', 'E', 'D', 'C', 'A', 'x'],
    'wz': ['y', 'E', 'B', 'w'],
}
time = [0 + m for m in range(50)]
departure = {'xy':0, 'zy':0, 'wy':0, 'yw':0, 'xz':0, 'xw':0, 'yz':0, 'yx':0, 'wx':0, 'wz':0, 'zx':0, 'zw':0}
arrival = {}
for f in flights:
    arrival[f] = departure[f]
    for s in paths[f]:
        arrival[f] += lengths[s]
    arrival[f] -= 1
depcap = {'x':1, 'y':1, 'z':1, 'w':1}
arrcap = {'x':1, 'y':1, 'z':1, 'w':1}
seccap = {'x':100, 'y':100, 'z':100, 'w':100, 'A':2, 'B':2, 'C':2, 'D':2, 'E':2, 'F':3}
ground = {'xy':1000, 'zy':1000, 'wy':1000, 'yw':1000, 'xz':1000, 'xw':1000, 'yz':1000, 'yx':1000, 'wx':1000, 'wz':1000, 'zx':1000, 'zw':1000}
air = {'xy':20, 'zy':20, 'wy':20, 'yw':20, 'xz':20, 'xw':20, 'yz':20, 'yx':20, 'wx':20, 'wz':20, 'zx':20, 'zw':20}

ghp = mip.Model(solver_name = "CBC")

x = [[[ghp.add_var(var_type=mip.BINARY, name=f'x_{f}_{t}_{j}') for j in range(len(paths[flights[f]]))] for t in range(len(time))] for f in range(len(flights))]

for k in airports:
    for t in range(1,len(time)):
        ghp += mip.xsum(
            (
                (x[f][t-1][0]-x[f][t][0])
                if paths[flights[f]][0] == k else 0
            )
            for f in range(len(flights))
        ) <= depcap[k]


for k in airports:
        for t in range(1,len(time)):
            ghp += mip.xsum(
                (
                    (x[f][t][len(paths[flights[f]])-1]-x[f][t-1][len(paths[flights[f]])-1])
                    if paths[flights[f]][-1] == k else 0
                )
                for f in range(len(flights))
            ) <= arrcap[k]

for j in sectors:
    for t in time:
        ghp += mip.xsum((x[f][t][paths[flights[f]].index(j)] if j in paths[flights[f]] else 0) for f in range(len(flights))) <= seccap[j]

for f in range(len(flights)):
    for t in time:
        ghp += mip.xsum(x[f][t][j] for j in range(len(paths[flights[f]]))) == 1

for f in range(len(flights)):
    for t in range(1,len(time)):
        for j in range(1,len(paths[flights[f]])):
            ghp += (x[f][t][j]-x[f][t-1][j]-x[f][t][j-1]+x[f][t-1][j-1]) >= -1

for f in range(len(flights)):
    for t in range(1,len(time)):
        for j in range(1,len(paths[flights[f]])):
            ghp += (2*x[f][t][j]-2*x[f][t-1][j]+x[f][t][j-1]-x[f][t-1][j-1]) <= 1

for f in range(len(flights)):
    for t in range(1,len(time)):
        for j in range(1,len(paths[flights[f]])):
            ghp += (x[f][t][j]-x[f][t-1][j]+2*x[f][t][j-1]-2*x[f][t-1][j-1]) >= -1

for f in range(len(flights)):
    for j in range(len(paths[flights[f]])):
        ghp += mip.xsum(x[f][t][j] for t in time) >= lengths[paths[flights[f]][j]]

for f in range(len(flights)):
    for t in range(departure[flights[f]]+1):
        ghp += x[f][t][0] == 1

for f in range(len(flights)):
        ghp += x[f][len(time)-1][len(paths[flights[f]])-1] == 1

for f in range(len(flights)):
    for t in range(len(time)-1):
        ghp += x[f][t+1][len(paths[flights[f]])-1] - x[f][t][len(paths[flights[f]])-1] >= 0

ghp.objective = mip.minimize(
    mip.xsum(
        (
            (ground[flights[f]]-air[flights[f]])
            *
            (
                mip.xsum(
                    ((t)*(x[f][t][0]-x[f][t+1][0]))for t in range(0,(len(time)-1))
                )
                -
                (departure[flights[f]])
            )
            +
            (air[flights[f]])
            *
            (
                mip.xsum(
                    (t*(x[f][t][len(paths[flights[f]])-1]-x[f][t-1][len(paths[flights[f]])-1]))for t in range(1,(len(time)))
                )
                -
                (arrival[flights[f]])
            )
        for f in range(len(flights)))
    )
)

solver_status = ghp.optimize()
objective_value= ghp.objective_value
print(f"Total cost: {objective_value}")

print("Time: ", end="")
for t in range(len(time)):
    if t<10:
        print("0"+str(t)+"  ", end="")
    else:
        print(str(t)+"  ", end ="")
print("\n")
for f in range(len(flights)):
    print("  "+flights[f]+":  ", end="")
    for t in range(len(time)):
        for j in range(len(paths[flights[f]])):
            var_value = ghp.var_by_name(f'x_{f}_{t}_{j}').x
            if var_value==1.0:
                print(paths[flights[f]][j] + "   ", end="")
        print("", end = "")
    print("\n")
            
