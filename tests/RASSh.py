''' N inputs A_1,A_2,...,A_N , N>=2 of fluid A with CFs, b_1,b_2,b_3,...,b_N.
    
    Pure buffer solution in stock. (N+1)>=3
	
    target CF = C_A
	
    integer d '''
import time
from z3 import *
import math as m
import os
#from graphviz import Digraph
from fractions import Fraction
counter =0

def getWeightedSum(dictionary,param_weight):
	'''This function takes the performance parameter in a dictionary and then 
	   outputs weighted sum of the parameters according to the param_weight list
		{'depth':depth,'n_m':n_m,'n_I':n_I,'n_b':n_b,'n_w':n_w}	
	'''
	value=0
	#parameters=dictionary.values()	
	#for i in range(len(parameters)):
	#	value+=parameters[i]*param_weight[i]
	value=param_weight[0]*dictionary['n_m']+param_weight[1]*dictionary['n_I']+param_weight[2]*dictionary['n_w']	
	return value

def optimal(param_weight,path):
	lst=[0]	
	for i in range(1,counter+1):
		f=open(os.path.join(path,'graph{}_parameters.txt'.format(i)),'r')
		dictionary=eval(f.readline())
		lst.append([getWeightedSum(dictionary,param_weight),dictionary])
		f.close()
	optimal_graphs=[]
	minimum=lst[1][0]
	for i in range(2,counter+1):
		if(minimum>lst[i][0]):
			minimum=lst[i][0]
	for i in range(1,counter+1):
		if minimum==lst[i][0]:
			if len(optimal_graphs)==0:
				optimal_graphs.append([lst[i],i])
			else:
				obj=optimal_graphs[0][0]
				if(obj[1]['n_m']<lst[i][1]['n_m'] or (obj[1]['n_m']==lst[i][1]['n_m'] and obj[1]['n_I']<lst[i][1]['n_I']) or (obj[1]['n_m']==lst[i][1]['n_m'] and obj[1]['n_I']==lst[i][1]['n_I'] and obj[1]['n_w']<lst[i][1]['n_w'])):
					optimal_graphs.pop(0)
					optimal_graphs.append([lst[i],i])  
				
	f=open(os.path.join(path,'optimal_graphs.txt'),'w')	 
	f.write(str(optimal_graphs))
	f.close()
	return optimal_graphs[0][0][1],minimum
		


def anyBaseToDecimal(lst,base):
	'''This method takes the input of a lst which contains the lst of the digits in the
	   number with radix as base. The MSB is in position lst[0]. Each of the digits in 
	   the representation in radix=base is an integer. The function returns the decimal
	   equivalent of this number in lst'''
	val=0
	n=len(lst)
	for i in range(0,n):
		val=val+(lst[i]*base**(n-i-1))	
	return val
	

def decimalToanyBase(val,base):
	'''This method takes as input a value in decimal and then it finds the radix=base
	   base equivalent of the value in lst with lst[0] containing the MSB, and each
	   of the digits in the lst is integer'''
	ac=''	
	while val!=0:
		ac=','+str(val%base)+ac
		val=val//base
	ac=ac[1:]
	if len(ac)!=0:
		lst=list(map(int,ac.split(sep=',')))
		return lst	
	else:
		return []		
	

def draw(X,D,d,N,B,m,path,M):
	depth=0	#height of the dilution tree
	n_m=0 	#total no. of (1:1) mixing steps
	n_I=0 	#total no. of input droplets (including buffer)
	n_b=0 	#total no. of buffer droplets
	n_w=0	#total no. of waste droplets
	global counter
	counter=counter+1
	k=1 #counts the number of nodes
	#pdf_with_path=os.path.join(path,'graph{}.gv'.format(counter))
	#g=Digraph('G',filename=pdf_with_path)
	#g.attr('node', shape='circle')
		
	lst=[]
	for j in range(d-1,-1,-1):
		for i in range(0,N):
			if X[i][j]!=0:
#				g.node(str(k),'A{}'.format(i+1)+':'+str(Fraction(B[i+1],M**d)))
#				g.node(str(k),shape='circle')
				lst.append([str(k),Fraction(B[i+1],M**d),X[i][j]])
				k=k+1
				if(depth==0):
					depth=j+1
				n_I+=X[i][j]
		if(D[j]!=0):
#			g.node(str(k),'D'+':'+str(Fraction(0,M**d)))
#			g.node(str(k),shape='circle')				
			lst.append([str(k),Fraction(0,M**d),D[j]])
			k=k+1
			#n_I+=D[j] (as n_I now calculates only the number of input reagents)
			n_b+=D[j]
		if(len(lst)==1):
			break
		q=[]
		while len(lst)!=0:
			#print('lst=',lst)
			val=0			
			index=0
			total=m
			while(total!=0):
				if lst[index][2]>total:
					val+=lst[index][1]*total
					total=0
				else:
					total=total-lst[index][2]
					val+=lst[index][1]*lst[index][2]
				index=index+1
			val=val/m
#			g.node(str(k),str(val))
#			g.node(str(k),style='filled',shape='circle',color='pink')
#			node=str(k)
			total=m
			for t in range(index):
				if lst[0][2]>total:
#					g.edge(lst[0][0],node,label=str(total))
					lst[0][2]-=total
					total=0
				else:
					total=total-lst[0][2]
#					g.edge(lst[0][0],node,label=str(lst[0][2]))	
					lst.pop(0)			
			q.append([str(k),val,1])
			k=k+1
			n_m+=1
		lst=q
#	g.node(lst[0][0],style='filled',shape='circle',color='green')
	numerator=lst[0][1].numerator
	denominator=lst[0][1].denominator
#	g.render(view=False)
	n_w=(m-1)*(n_m-1)
	f1=open(os.path.join(path,'graph{}_parameters.txt'.format(counter)),'w+')
	dictionary={'depth':depth,'n_m':n_m,'n_I':n_I,'n_b':n_b,'n_w':n_w,'exact':lst[0][1],'num':numerator,'deno':denominator}
	#'C_t':(numerator/(denominator/m**depth)),
	f1.write(str(dictionary)) # use eval(f.readline()) to read the dictionary from file	
	f1.close()
					
	

def display(a):
	for i in a:
		for j in i:
			print(j,end='\t')
		print()

def process(X,d,N,B,m,path,M):
		
	count=0
	for i in range(len(X)):
		count=count+anyBaseToDecimal(X[i],m)
	count=m**d-count
	beta=decimalToanyBase(count,m)
	beta=[0]*(d-len(beta))+beta		
	#print()
	#print(beta)
	#print()
	draw(X,beta,d,N,B,m,path,M)
'''
1st column =  M (mixer size)
2nd column = d (depth of the graph)
3rd column = R (number of  arbitrary stock solutions)
4th to 6th column = A_1 A_2 A_3 (Set of arbitrary stock solutions)
7th column = C_t (Target concentration)
8th column = M^d (denominator of target concentration) '''
def generate(lst,path):

	start=time.time() # get the start time
	global counter
	counter=0
	print("THIS IS THE MEDA-HEURISTIC ALGORITHM FOR N>=2")
	M=lst[0]#number of mixing models
	d=lst[1] #int(input("ENTER THE VALUE OF d : "))
	N=lst[2]#int(input("ENTER THE VALUE OF N : "))
	#print("ENTER THE CFs OF EACH OF THE INPUT FLUIDS...")
	param_weight=[Fraction(6,10),Fraction(3,10),Fraction(1,10)]
	#{'depth':depth,'n_m':n_m,'n_I':n_I,'n_b':n_b,'n_w':n_w}

	B=[0]

	for i in range(N):
		val=lst[3+i] #float(input("ENTER THE CF OF FLUID {} ".format(i+1)))
		B.append(val)

	B.append(0) # for the buffer

	T=lst[-2]
	#T_rational=Q(T,1)
	#C_A=(f1.readline()) #float(input("ENTER THE TARGET CF : "))



	deno=lst[-1]

	#B=[0]

	#for i in range(1,N+2):

	#T=m.floor(C_A*2**d)
	for m in range(M,1,-1):
		print(" m = ",m)
		W=[[Q(B[i],m**j) for j in range(0,d+1)] for i in range(0,N+2)]

		wt=[Q(1,m**j) for j in range(0,d+1)]

		X=[[Int('X{}_{}'.format(i,j)) for j in range(0,d+1)]for i in range(0,N+2)]

		V=[]

		T1=[]

		s=Solver()

		for i in range(1,N+1):
			for j in range(1,d+1):
				T1.append(X[i][j]*W[i][j])
				V.append(X[i][j]*wt[j])
				s.add(And(0<=X[i][j],X[i][j]<=m-1))
			
		s.add(And((T-0.5)<Sum(T1),Sum(T1)<(T+0.5)))
		s.add(Sum(V)<=1)

		while s.check()==sat:
			if time.time()>start+120:
				break
			model=s.model()
			subset=[]
			cond=[]
			for i in range(1,N+1):
				lst=[]
				for j in range(1,d+1):
					lst.append(int(str(model[X[i][j]])))
					cond.append(model[X[i][j]]==X[i][j])
				subset.append(lst)
			s.add(Not(And(cond)))
			#print()
			#display(subset)
			#print()
			process(subset,d,N,B,m,path,M)
			#print()
	'''dictionary={'depth':depth,'n_m':n_m,'n_I':n_I,'n_w':n_w,'C_t':(numerator/(denominator/m**depth)),'exact':lst[0][1]}'''
	dictionary,minimum=optimal(param_weight,path)
	depth=dictionary['depth']
	n_m=dictionary['n_m']
	n_I=dictionary['n_I']
	n_w=dictionary['n_w']
	n_b=dictionary['n_b']
	numerator=dictionary['num']
	denominator=dictionary['deno']	
	C_t=numerator/(denominator/deno)
	exact=dictionary['exact']
	end =time.time()
	f=open(os.path.join(path,'time.txt'),'w')
	f.write(str(end-start)+' s')
	f.close()
	r_t='%.4f'%(end-start)
	n_sol=counter
	delta='%.4f'%(abs(C_t-T))#abs(float(exact-Fraction(T,M**d)))
	C_t='%.4f'%C_t
	minimum='%.4f'%minimum
	string=str(depth)+'\t'+str(C_t)+'\t'+str(deno)+'\t'+str(n_m)+'\t'+str(n_I)+'\t'+str(n_b)+'\t'+str(n_w)+'\t'+str(r_t)+'\t'+str(n_sol)+'\t'+str(delta)+'\t'+str(minimum)
	return string #,[n_m,n_I,n_b,n_w,f]
	'''d Ct' M^d n_m n_I n_b n_w r_t n_sol delta f'''

def readInput():
	f=open('input.txt','r')
	filename=f.readline()
	while len(filename)!=0:
		filename=filename[:-1]	#to remove the extra \n read at the end
		f1=open(filename,'r')
		f2=open(filename[:-4]+'_out.txt','w')		
		os.makedirs(filename[:-4])
		line=f1.readline()
		i=1
		accumulator=[0,0,0,0]
		while len(line)!=0:		
			lst=list(map(int,line.split()))
			path=os.path.join(filename[:-4],str(i)+'.'+str(lst[-2])+'_'+str(lst[-1]))
			os.makedirs(path)
			string=generate(lst,path)
			f2.write(str(string)+'\n')
			line=f1.readline()
			i=i+1
			#for j in range(0,5):
			#	accumulator[j]+=lst[j]
		#for j in range(0,5):
		#	accumulator[j]/=(i-1)
		#f3=open(os.path.join(filename[:-4],'average.txt'),'w')
		#f3.write(str(accumulator))
		#f3.close()			
		f1.close()
		f2.close()
		filename=f.readline()
	f.close()

readInput()


