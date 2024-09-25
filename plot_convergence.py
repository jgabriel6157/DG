import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

x = [15,31,63,127]

#Order 1
t1 = [460,873,1706,3440.5,6790]
y1 = [0.0708,0.00303,0.00176,0.00171555]
# x = [10]

#Order 2
t2 = [98,163.5,292.5,563.7,1063.86]
y2 = [0.07108,0.0062623,0.005621,0.005577]
# x = [4,8,16,32,64]

#Order 4
t4 = [241,439.8,851,1690,3343]
y4 = [0.00429,0.000126,5.237e-6,9.7589e-7]
# x = [4,8,16,32,64]

#Order 8
t8 = [728,1417,2847,5662,11152]
y8 = [7.887e-6,1.43e-6,6.965e-7,4.482e-7]
# x = [4,8,16,32,64]

plt.plot(x,y1,marker='s')
# plt.plot(x,y2,marker='s')
# plt.plot(x,y4,marker='s')
# plt.plot(x,y8,marker='s')

# plt.title("Order = 2")
plt.xlabel("Nz")
# plt.xlabel("Nz")
plt.ylabel("L2 Error Norm")
# plt.legend(["N=1","N=2","N=4","N=8"])
plt.ylim((0.001,0.1))
plt.yscale("log")
plt.show()