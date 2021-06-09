import matplotlib.pyplot as plt
def main():
	i = 0
	f = open('C:/Users/User/source/repos/CHW_9/CHW_9/linear.txt','r')
	l = []
	for line in f:
    		l.append(list(map(float, line.split())))
    		i+=1
	plt.plot(l[0],l[1],label = 'ballistic')
	plt.plot(l[0],l[2],label = 'finite_substrations')
	plt.legend()
	plt.show()
main()