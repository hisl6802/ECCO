import logging
import time
import os
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import tkinter as tk
import multiprocessing
import webbrowser
import GuiBackground as GB
from GUIUtils import GUIUtils as GU
import config


class JuneLabClusteringGUI(ttk.Frame):
	def __init__(self, master=None):
		super().__init__(master)
		self.grid(column=0, row=0, sticky=(N, W, E, S))
		self.rowconfigure(0, weight =2)
		self.columnconfigure(1, weight =2)
		self.columnconfigure(2, weight =2)
		self.rowconfigure(1, weight = 2)
		self.rowconfigure(2, weight = 2)
		self.rowconfigure(3, weight = 2)
		self.rowconfigure(4, weight = 2)
		self.columnconfigure(3, weight=2)
	
		self.startUpPage()

	def startUpPage(self):
		'''
		Working on a new starting page to ensure that the user selects the appropriate number of threads, can input there name etc.
		'''
		self.style = ttk.Style()
		self.style.configure("RW.TLabel", foreground="#f03d33",font=("TkHeadingFont",30))
		self.style.configure("RW.TButton", padding=15, borderwidth=15, foreground="black", background="#000000",font=("Arial",14))

		numThreads = int(multiprocessing.cpu_count())
		#set up the start up page.
		self.JuneLab = ttk.Label(self, text="GUI Set-Up",style="RW.TLabel").grid(column=0,row=0,columnspan=4)
		self.NameLab = ttk.Label(self, text="Please input your name or a Project name:",font=("TkHeadingFont",16)).grid(column=2,row=1,sticky=(N,S,E,W),pady=10)
		self.name = tk.StringVar()
		self.threads = tk.StringVar()
		self.entryName = ttk.Entry(self,textvariable=self.name).grid(column=2,row=2,sticky=(N,S,E,W),pady=3,padx=5)

		self.ThreadsLab = ttk.Label(self, text="Number of threads:",font=("TkHeadingFont",16)).grid(column=2,row=3,sticky=(N,S,E,W),pady=1)
		self.entryThreads = ttk.Entry(self,textvariable=self.threads).grid(column=2,row=4,sticky=(N,S,E,W), pady=1,padx = 5)
		self.ThreadsULab = ttk.Label(self, text="You have "+str(numThreads)+ " available. (Using half or less is recommended.)",font=("TkHeadingFont",16)).grid(column=2,row=5,sticky=(N,S,E,W),pady=2)
		self.getStarted = ttk.Button(self,text="Get Started!",command=self.create_widgets).grid(column=2, row=6,sticky=(N,S,E,W),columnspan=1)

	def create_widgets(self):
		'''
		'''

		#get the project name
		name = self.name.get()
		numThreads = self.threads.get()
		if int(numThreads) <= multiprocessing.cpu_count():
			config.numThreads = int(numThreads)
		else:
			config.numThreads = 2
			messagebox.showinfo(title="Number of Threads", message="You have been assigned 2 threads, since an invalid number of input threads was detected.")

		config.name = name

		#go to the appropriate directory and create log file string
		log_time = time.strftime("%a_%b_%d_%Y_%H_%M_%S")
		log_file = config.name + '_' + str(log_time) + '.log' 

		logging.basicConfig(filename=log_file,format='%(asctime)s %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO)
		logging.info('Started GUI')

		objects = self.grid_slaves()
		for i in objects:
			i.destroy()

		self.style = ttk.Style()
		self.style.configure("RW.TLabel", foreground="#f03d33",font=("TkHeadingFont",30))
		self.style.configure("RW.TButton", padding=15, borderwidth=15, foreground="black", background="#000000",font=("Arial",14))
		self.JuneLab = ttk.Label(self, text="ECCO GUI",style="RW.TLabel").grid(column=0,row=0,columnspan=4)
		self.clust = ttk.Button(self,text="Create Clustergram",style="RW.TButton",command=self.createClustergram).grid(column=1,row=1, sticky=(N,S,E,W))
		#Create a button to allow the user to create a medians file for better clustering results. 
		self.med = ttk.Button(self, text="Group Medians", style="RW.TButton", command=self.medians).grid(column=1, row=3, sticky =(N,S,E,W))
		#Create a button to allow the user to compare the four most common linkage functions.
		self.link = ttk.Button(self, text="Compare Linkage Functions", style="RW.TButton", command=self.linkages).grid(column=2, row=1,sticky =(N,S,E,W))
		#Create a button to allow the user to create the peaks to pathways files needed to analyze the peaks to pathways in Mummichog
		self.peak = ttk.Button(self, text="Peaks to Pathways", style="RW.TButton", command=self.P2P).grid(column=2, row=2, sticky =(N,S,E,W))
		#Create a button to allow the user to check the integrity of the data downloaded from Metaboanalysts Volcano plot results. 
		self.integrity = ttk.Button(self, text="Data Integrity", style="RW.TButton", command=self.integrity).grid(column=2, row=3, sticky =(N,S,E,W))
		#Create a button to allow the user to do an Ensemble clustering on the data.
		self.ensemble = ttk.Button(self,text="Ensemble Clustering", style="RW.TButton", command=self.ensemble).grid(column=1,row=2,sticky=(N,S,E,W))
		#Create a button to allow the user to create a minimum spanning tree on data
		self.mst = ttk.Button(self,text='Cluster Optimization', style="RW.TButton", command=self.mst).grid(column=3,row=1,sticky=(N,S,E,W))
		#Create a button for the generation of a report
		self.generate = ttk.Button(self,text='Selected Clusters Figure', style="RW.TButton", command=self.genSelClustFig).grid(column=3,row=3,sticky=(N,S,E,W))
		#Create a button for the selection of clusters
		self.selection = ttk.Button(self, text="Cluster Selection",style = "RW.TButton",command=self.clusterSelection).grid(column=3, row=2,sticky=(N,S,E,W))
		# pad each widget with 5 pixels on each side to ensure that the buttons do not stay together. 
		for child in self.winfo_children(): child.grid_configure(padx=5, pady=5)

	def home(self):
		#start listing out the global variables that need to be removed for each function upon returning home
		if 'ensemble' in globals():
			try:
				del(globals()['ensemble'])
			except:
				messagebox.showerror(title='Error', message="GUI didn't properly reset. Restart recommended, IF you want to run another ensemble.")
				logging.error(': Ensemble environment did not delete, this may cause errors in the computation of wanted ensemble clustergram!')
			
			if 'standard' in globals():
				try:
					del(globals()['standard'])
				except:
					messagebox.showerror(title='Error', message="GUI didn't properly reset. Restart recommended, IF you want to run another ensemble.")
					logging.error(': Standard ensemble was selected but global variable didn''t delete this may cause issues with other ensemble clusterings!')

		#remove old objects and put the home objects back on grid
		objects = self.grid_slaves()
		for i in objects:
			i.grid_remove()
		widgets = self.winfo_children()

		widgetDict = {}
		for i in range(10):
			#create a dictionary of the widgets from home window
			widgetDict[i] = widgets[i]

		widgetDict[0].grid(column=0, row=0,columnspan=4)
		widgetDict[1].grid(column=1, row=1, sticky=(N,S,E,W))
		widgetDict[2].grid(column=1, row=3, sticky =(N,S,E,W))
		widgetDict[3].grid(column=2, row=1,sticky =(N,S,E,W))
		widgetDict[4].grid(column=2, row=2, sticky =(N,S,E,W))
		widgetDict[5].grid(column=2, row=3, sticky =(N,S,E,W))
		widgetDict[6].grid(column=1, row=2, sticky =(N,S,E,W))
		widgetDict[7].grid(column=3, row=1,sticky=(N,S,E,W))
		widgetDict[8].grid(column=3, row=3,sticky=(N,S,E,W))
		widgetDict[9].grid(column=3,row=2,sticky=(N,S,E,W))


		count = -1
		for child in self.winfo_children():
			#add padding to the current widgets
			count += 1
			if count < 10:
				child.grid_configure(padx=5,pady=5)


	def createClustergram(self):
		def linkageOutput(*args):
			#grab the current selection of the list
			global selection
			selection = distListBox.curselection()
			curLink = linkageList[selection[0]]
			if curLink == 'ward':
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				sampleListBox.insert(0,distList[0])
				
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			else:
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				for i in range(len(distList)):
					sampleListBox.insert(i,distList[i])

				#bind the output back to the GUI.
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			return selection

		def colorMap(*args):
			global selection1
			selection1 = sampleListBox.curselection()
			#create a list of the color map options
			lenList = len(colorListBox.get(0,tk.END))
			if lenList > 0:
				colorListBox.delete(0,lenList-1)

			for i in range(len(colorList)):
				colorListBox.insert(i,colorList[i])

			#bind the output back to the GUI
			self.colorListBox = colorListBox
			self.colorListBox.grid(column=3,row=2,columnspan=1)
			return selection1

		def dataTransform(*args):
			#dataTransform
			global selection2
			selection2 = colorListBox.curselection()

			#put the data transform options into the list
			lenList = len(transformListBox.get(0,tk.END))
			if lenList > 0:
				transformListBox.delete(0,lenList-1)

			for i in range(len(transformList)):
				transformListBox.insert(i,transformList[i])

			self.transformListBox = transformListBox
			self.transformListBox.grid(column=1,row=4,columnspan=1)
			return selection2

		def dataScale(*args):
			#dataScaling
			global selection3
			selection3 = transformListBox.curselection()
			
			#put the data scaling optoins into the list
			lenList = len(scaleListBox.get(0,tk.END))
			if lenList > 0:
				scaleListBox.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleListBox.insert(i,scaleList[i])

			self.scaleListBox = scaleListBox
			self.scaleListBox.grid(column=2,row=4,columnspan=1)
			return selection3

		def dataNorm(*args):
			#Does the user want to normalize to a column? 
			global selection4
			selection4 = scaleListBox.curselection()

			#put the options into the list
			#put the data scaling optoins into the list
			lenList = len(normListBox.get(0,tk.END))
			if lenList > 0:
				normListBox.delete(0,lenList-1)

			for i in range(len(normList)):
				normListBox.insert(i,normList[i])

			self.normListBox = normListBox
			self.normListBox.grid(column=3,row=4,columnspan=1)


		def submit(*args):
			#submit the selections to the function output
			selection5 = normListBox.curselection()

			dist = distList[selection1[0]]
			link = linkageList[selection[0]]
			color = colorList[selection2[0]]
			transform = transformList[selection3[0]]
			scale = scaleList[selection4[0]]
			norm = normList[selection5[0]]

			if norm == 'Normalize':
				scale = 'NormStand'
				groupOrd = inputGroupOrder.get()
				groupOrd = groupOrd.split(',')
				GU.createClustergram(1,link,dist,color,colOrder=groupOrd, transform=transform,scale=scale)
			else:
				groupOrd = inputGroupOrder.get()
				groupOrd = groupOrd.split(',')
				norm = 0
				if len(groupOrd) > 1:
					norm =2

				GU.createClustergram(norm,link,dist,color,colOrder=groupOrd, transform=transform,scale=scale)

		def cmapO(*args):
			#send users to webpage of 
			webbrowser.open('https://matplotlib.org/stable/tutorials/colors/colormaps.html')

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Clustergram Input",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.Linkage = ttk.Label(self, text="Linkage",font=("TkHeadingFont",12)).grid(column=1,row=1,sticky=(N))
		self.Distance = ttk.Label(self, text="Distance Measure",font=("TkHeadingFont",12)).grid(column=2,row=1,sticky=(N))
		self.Color = ttk.Label(self,text="Color-Map", font=("TkHeadingFont",12)).grid(column=3,row=1,sticky=(N))
		self.Transform = ttk.Label(self,text="Transform", font=("TkHeadingFont",12)).grid(column=1,row=3,sticky=(N))
		self.Scale = ttk.Label(self,text="Scale",font=('TkHeadingFont',12)).grid(column=2,row=3,sticky=(N))
		self.homepage = ttk.Button(self,text="Return to Home",command=self.home).grid(column=2, row=8,sticky=(N),columnspan=1)
		self.submit = ttk.Button(self,text="Submit", command=submit).grid(column=2, row=7,sticky=(N),columnspan=1)
		self.cmapW = ttk.Button(self,text="ColorMap Options", command=cmapO).grid(column=3,row=3,sticky=(N),columnspan=1)
		self.gLab = ttk.Label(self,text="Group Order, normilization column first",font=('TkHeadingFont',12)).grid(column=1,row=5,columnspan=3)
		inputGroupOrder = tk.StringVar()
		self.inputGroups = ttk.Entry(self,textvariable=inputGroupOrder).grid(column=2,row=6,sticky=(N))
		distListBox = Listbox(self,height=8)
		sampleListBox = Listbox(self,height=8)
		colorListBox = Listbox(self,height=8)
		transformListBox = Listbox(self, height=8)
		scaleListBox = Listbox(self, height=8)
		normListBox = Listbox(self,height=8)
		
		#Create the lists of available options for selection 
		linkageList = config.linkageList
		distList = config.distList
		colorList = config.colorList
		transformList = config.transformList 
		scaleList = config.scaleList 
		normList = config.normList

		
		linkNames = StringVar(value=linkageList)
		distNames = StringVar(value=distList)
		colorNames = StringVar(value=colorList)
		transformNames = StringVar(value=transformList)
		scaleNames = StringVar(value=scaleList)
		normNames = StringVar(value=normList)

		
		#input the linkage function values into the box
		for i in range(len(linkageList)):
			distListBox.insert(i,linkageList[i])

		distListBox.bind('<Double-1>',linkageOutput)
		sampleListBox.bind('<Double-1>',colorMap)
		colorListBox.bind('<Double-1>',dataTransform)
		transformListBox.bind('<Double-1>', dataScale)
		scaleListBox.bind('<Double-1>',dataNorm)
		self.distListBox = distListBox
		self.distListBox.grid(column=1,row=2,columnspan=1)
		self.sampleListBox = sampleListBox
		self.sampleListBox.grid(column=2,row=2,columnspan=1)
		self.colorListBox = colorListBox
		self.colorListBox.grid(column=3,row=2,columnspan=1)
		self.transformListBox = transformListBox
		self.transformListBox.grid(column=1,row=4,columnspan=1)
		self.scaleListBox = scaleListBox
		self.scaleListBox.grid(column=2,row=4,columnspan=1)
		self.normListBox = normListBox
		self.normListBox.grid(column=3,row=4,columnspan=1)

	def clusterSelection(self):
		def linkageOutput(*args):
			#grab the current selection of the list
			global selection
			selection = distListBox.curselection()
			curLink = linkageList[selection[0]]
			if curLink == 'ward':
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				sampleListBox.insert(0,distList[0])
				
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			else:
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				for i in range(len(distList)):
					sampleListBox.insert(i,distList[i])

				#bind the output back to the GUI.
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			return selection

		def colorMap(*args):
			global selection1
			selection1 = sampleListBox.curselection()
			#create a list of the color map options
			lenList = len(colorListBox.get(0,tk.END))
			if lenList > 0:
				colorListBox.delete(0,lenList-1)

			for i in range(len(colorList)):
				colorListBox.insert(i,colorList[i])

			#bind the output back to the GUI
			self.colorListBox = colorListBox
			self.colorListBox.grid(column=3,row=2,columnspan=1)
			return selection1

		def dataTransform(*args):
			#dataTransform
			global selection2
			selection2 = colorListBox.curselection()

			#put the data transform options into the list
			lenList = len(transformListBox.get(0,tk.END))
			if lenList > 0:
				transformListBox.delete(0,lenList-1)

			for i in range(len(transformList)):
				transformListBox.insert(i,transformList[i])

			self.transformListBox = transformListBox
			self.transformListBox.grid(column=1,row=4,columnspan=1)
			return selection2

		def dataScale(*args):
			#dataScaling
			global selection3
			selection3 = transformListBox.curselection()
			#put the data scaling optoins into the list
			lenList = len(scaleListBox.get(0,tk.END))
			if lenList > 0:
				scaleListBox.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleListBox.insert(i,scaleList[i])

			self.scaleListBox = scaleListBox
			self.scaleListBox.grid(column=2,row=4,columnspan=1)
			
			return selection3

		def dataNorm(*args):
			#Does the user want to normalize to a column? 
			global selection4
			selection4 = scaleListBox.curselection()

			#put the options into the list
			#put the data scaling optoins into the list
			lenList = len(normListBox.get(0,tk.END))
			if lenList > 0:
				normListBox.delete(0,lenList-1)

			for i in range(len(normList)):
				normListBox.insert(i,normList[i])

			self.normListBox = normListBox
			self.normListBox.grid(column=3,row=4,columnspan=1)
			self.submit.grid(column=3, row=8,sticky=(N),columnspan=1)


		def submit(*args):
			#submit the function output to the
			norm = normListBox.curselection()
			dist = distList[selection1[0]]
			link = linkageList[selection[0]]
			color = colorList[selection2[0]]
			transform = transformList[selection3[0]]
			scale = scaleList[selection4[0]]
			norm = normList[norm[0]]

			#set the config colorNum to zero
			config.colorNum = 0
			if norm == 'Normalize':
				groupOrd = inputGroupOrderCS.get()
				groupOrd = groupOrd.split(',')
				scale ='NormStand'
				GU.selectClusters(link,dist,1,colOrder=groupOrd,transform=transform, scale=scale,cmap=color)
			else:
				GU.selectClusters(link,dist,transform=transform, scale=scale,cmap=color)



		def cmapO(*args):
			#send users to webpage of 
			webbrowser.open('https://matplotlib.org/stable/tutorials/colors/colormaps.html')

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()


		
		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Clustergram Input",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.Linkage = ttk.Label(self, text="Linkage",font=("TkHeadingFont",12)).grid(column=1,row=1,sticky=(N))
		self.Distance = ttk.Label(self, text="Distance Measure",font=("TkHeadingFont",12)).grid(column=2,row=1,sticky=(N))
		self.Color = ttk.Label(self,text="Color-Map", font=("TkHeadingFont",12)).grid(column=3,row=1,sticky=(N))
		self.Transform = ttk.Label(self,text="Transform", font=("TkHeadingFont",12)).grid(column=1,row=3,sticky=(N))
		self.Scale = ttk.Label(self,text="Scale",font=('TkHeadingFont',12)).grid(column=2,row=3,sticky=(N))
		self.homepage = ttk.Button(self,text="Return to Home",command=self.home).grid(column=2, row=8,sticky=(N),columnspan=1)
		self.submit = ttk.Button(self,text="Submit", command=submit)
		self.cmapW = ttk.Button(self,text="ColorMap Options", command=cmapO).grid(column=1,row=8,sticky=(N),columnspan=1)
		self.groupsLAb = ttk.Label(self,text="Group Order, normilization column first",font=('TkHeadingFont',12)).grid(column=2,row=6,sticky=(N))
		inputGroupOrderCS = tk.StringVar()
		self.inputGroups = ttk.Entry(self,textvariable=inputGroupOrderCS).grid(column=2,row=7,sticky=(N))
		distListBox = Listbox(self,height=8)
		sampleListBox = Listbox(self,height=8)
		colorListBox = Listbox(self,height=8)
		transformListBox = Listbox(self, height=8)
		scaleListBox = Listbox(self, height=8)
		normListBox = Listbox(self, height=8)
		
		#Create the lists of available options for selection 
		linkageList = ('single','ward','complete','average')
		distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
		colorList = ('viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','twilight', 'twilight_shifted', 'hsv')
		transformList = ('None','Log transformation', 'Square root transformation', 'Cube root transformation')
		scaleList = ('None', 'Mean centering', 'Auto Scaling', 'Pareto Scaling', 'Range Scaling')
		normList = config.normList

		
		linkNames = StringVar(value=linkageList)
		distNames = StringVar(value=distList)
		colorNames = StringVar(value=colorList)
		transformNames = StringVar(value=transformList)
		scaleNames = StringVar(value=scaleList)
		normNames  = StringVar(value=normListBox)

		
		#input the linkage function values into the box
		for i in range(len(linkageList)):
			distListBox.insert(i,linkageList[i])

		distListBox.bind('<Double-1>',linkageOutput)
		sampleListBox.bind('<Double-1>',colorMap)
		colorListBox.bind('<Double-1>',dataTransform)
		transformListBox.bind('<Double-1>', dataScale)
		scaleListBox.bind('<Double-1>', dataNorm)
		self.distListBox = distListBox
		self.distListBox.grid(column=1,row=2,columnspan=1)
		self.sampleListBox = sampleListBox
		self.sampleListBox.grid(column=2,row=2,columnspan=1)
		self.colorListBox = colorListBox
		self.colorListBox.grid(column=3,row=2,columnspan=1)
		self.transformListBox = transformListBox
		self.transformListBox.grid(column=1,row=4,columnspan=1)
		self.scaleListBox = scaleListBox
		self.scaleListBox.grid(column=2,row=4,columnspan=1)
		self.normListBox =normListBox
		self.normListBox.grid(column=3,row=4,columnspan=1)

	def medians(self):
		global rmZeros
		#function to 
		def groupMedians(*args):
			rmZeros = var1.get()
			GU.groupMedians(rmZeros=rmZeros)
		#get rid of the objects in GUI window.
		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the group medians.
		self.GroupLab = ttk.Label(self, text="Group Medians", font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=2)
		self.homeGroup = ttk.Button(self,text="Return to Home",command=self.home).grid(column=1,row=4,sticky=(N),columnspan=2)
		var1 = IntVar()
		self.zeroRemove = ttk.Checkbutton(self,text="Remove Zeros?",variable=var1).grid(column=1,row=2,sticky=(N),columnspan=2)

		self.GM = ttk.Button(self,text="Select file",command=groupMedians).grid(column=1,row=3,sticky=(N),columnspan=2)


	def linkages(self):
		def distFunc(*args):
			#make a global distance variable
			global distanceMet
			distanceMet = self.dist.get()
			
			#given the selection of a distance measure how many comparisons are possible. 
			if distanceMet == 'euclidean':
				#give the full list to the second combobox
				values = [1,2,3,4]
				num_comps = StringVar()
				self.numComps = ttk.Combobox(self,values=values,textvariable=num_comps)
			else:
				#give only three values to the second combobox
				values = [1,2,3]
				num_comps = StringVar()
				self.numComps = ttk.Combobox(self,values=values,textvariable=num_comps)
			
			#place the combobox on the GUI
			self.numComps.bind('<<ComboboxSelected>>', numCompsFunc)
			self.numComps.grid(column=2,row=2)
		
		def numCompsFunc(*args):
			#make a numberComps global variable
			global numberComps
			numberComps = self.numComps.get()

			#create a list of linkage options
			linkOpts = ['ward-single','ward-complete','ward-average','single-complete','single-average','complete-average',\
						'ward-single-complete','ward-single-average','ward-complete-average','single-complete-average',\
						'ward-single-complete-average']
			
			if numberComps == '1':
				#give the linkage funtions combobox a list of the linkage functions. 
				linkages = ['ward','single','complete','average']

				#check for the euclidean measure again.
				if distanceMet == 'euclidean':
					#give the values of linkages[1:3]
					value = linkages
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				else:
					value = linkages[1:]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values = value,textvariable=linkage)

				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

			elif numberComps == '2':
				#give the linkage functions combobox a subset based upon distance measure. 
				if distanceMet == 'euclidean':
					#give the list values from 0 to 5
					value = linkOpts[0:6]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)

				else:
					#give the list values from 3 to 5
					value = linkOpts[3:6]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				
				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

			elif numberComps == '3':
				#give the linkage functions combobox a subset based upon distance measure. 
				if distanceMet == 'euclidean':
					#give the list values from 0 to 5
					value = linkOpts[6:10]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)

				else:
					#give the list values from 3 to 5
					value = linkOpts[9]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				
				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

			elif numberComps == '4':
				#give the list values from 3 to 5
				value = linkOpts[10]
				linkage = StringVar()
				self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				
				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

		def linkageComp(*args):
			#base linkage list for 4 comparisons.
			global linkList
			linkageList = ['ward','single','complete','average']
			
			#create an empty list for linkage functions
			linkList = []

			#count the number of dashes in the string from combobox
			selection = self.linkage.get()
			locs = []
			for i in range(len(selection)):
				#find the dash locations
				if selection[i] == '-':
					locs.append(i)

			if len(locs) > 0:
				firstLetter = 0
				for i in range(len(selection)):
					#check vthe current string value for '-'
					if selection[i] == '-':
						lastLetter = i
						#get current linkage
						curLink = selection[firstLetter:lastLetter]
						linkList.append(curLink)
						firstLetter = i+1
					elif i == len(selection)-1:
						lastLetter = i+1
						curLink = selection[firstLetter:lastLetter]
						linkList.append(curLink)

			else:
				#append the selection to the linkage list
				linkList.append(selection)
			lenList = len(transformListBox.get(0,tk.END))
			if lenList > 0:
				transformListBox.delete(0,lenList-1)
			#send the parameters for linkage comparison 
			for i in range(len(transformList)):
				transformListBox.insert(i,transformList[i])	

		def dataScale(*args):
			global dataTrans
			dataTrans = transformListBox.curselection()
			dataTrans = transformList[dataTrans[0]]
			lenList = len(scaleListBox.get(0,tk.END))
			if lenList > 0:
				scaleListBox.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleListBox.insert(i,scaleList[i])
			
		def submit(*args):
			#get current selection of the scaling list box
			dataScale = scaleListBox.curselection()
			dataScale = scaleList[dataScale[0]]

			file = filedialog.askopenfilename()
			GU.linkageComparison(file, numberComps,linkList,distanceMet,dataTrans,dataScale)


		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create a list of values from 1 to 4
		numLinkComps = [1,2,3,4]

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Linkage Comparison",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.numCompsLab = ttk.Label(self, text="Number of comparisons",font=("TkHeadingFont",12)).grid(column=2,row=1)
		self.distLab = ttk.Label(self, text="Distance measure",font=("TkHeadingFont",12)).grid(column=1,row=1)
		self.linkLab = ttk.Label(self,text="Linkage functions",font=("TkHeadingFont",12)).grid(column=3,row=1)
		self.Transform = ttk.Label(self,text="Transform",font=("TkHeadingFont",12)).grid(column=1,row=4)
		self.Scale = ttk.Label(self,text="Scale",font=("TkHeadingFont",12)).grid(column=2,row=4)
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=1,row=7, sticky=(N),columnspan=3)
		self.sumbitIt = ttk.Button(self,text="Submit",command=submit).grid(column=1,row=6, sticky=(N),columnspan=3)
		
		transformListBox = Listbox(self, height=8)
		scaleListBox = Listbox(self, height=8)

		linkages = StringVar()
		#create the distance measure combobox first, then update the GUI as the user selects the distance, measure than number of linkage comps. 
		distances = StringVar()
		distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
		transformList = ('None','Log transformation', 'Square root transformation', 'Cube root transformation')
		scaleList = ('None', 'Mean centering', 'Auto Scaling', 'Pareto Scaling', 'Range Scaling')
		self.dist = ttk.Combobox(self,values = distList,textvariable=distances)


		
		#bind the combobox for distance measures to the selection of distance measure. 
		self.dist.bind('<<ComboboxSelected>>', distFunc)
		self.dist.grid(column=1,row=2)
		transformListBox.bind('<Double-1>', dataScale)
		scaleNames = StringVar(value=scaleList)
		self.transformListBox = transformListBox
		self.transformListBox.grid(column=1,row=5,columnspan=1)
		self.scaleListBox = scaleListBox
		self.scaleListBox.grid(column=2,row=5,columnspan=1)

	def P2P(self):
		GU.peaksToPathways()

	def integrity(self):
		filename = filedialog.askopenfilename()
		GU.dataIntegrity()
	
	def mstF(self):
		numClust =GU.MST(func='ensemble')
		numClust = int(numClust)
		GU.ensembleClustering(optNum=numClust)

	def ensemble(self):
		global ensemble
		ensemble = 1

		def firstNext(*args):
			global ensemLinkList 
			singleBoxOut = singleBox.get()
			completeBoxOut = completeBox.get()
			averageBoxOut = averageBox.get()
			wardBoxOut = wardBox.get()

			linkList = [singleBoxOut,completeBoxOut,averageBoxOut,wardBoxOut]

			ensemLinkList = []
			for i in range(len(linkList)):
				if len(linkList[i]) > 0:
					ensemLinkList.append(linkList[i])


			## Selection error handling---------------------------------------------------------------------------------------------------------------------------------------------------------
			if len(ensemLinkList) == 0:
				messagebox.showerror(title="Seletion Error", message="You must select at least one linkage function!")
				return

			elif len(ensemLinkList) == 1:

				if 'ward' in ensemLinkList:
					messagebox.showerror(title="Selection Error", message="The ward linkage function only works with the Euclidean distance function. Thus, you cannot run an ensemble. Please select at least one other linkage function if you would like to use ward.")
					return
			##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			objects = self.grid_slaves()
			for i in objects:
				i.grid_forget()
			global euclideanBox
			global seuclideanBox
			global sqeuclideanBox
			global cosineBox
			global chebyshevBox
			global correlationBox
			global canberraBox
			global braycurtisBox
			global minkowskiBox
			global cityblockBox

			euclideanBox = tk.StringVar()
			seuclideanBox = tk.StringVar()
			sqeuclideanBox = tk.StringVar()
			cosineBox = tk.StringVar()
			chebyshevBox = tk.StringVar()
			correlationBox = tk.StringVar()
			canberraBox = tk.StringVar()
			braycurtisBox = tk.StringVar()
			minkowskiBox = tk.StringVar()
			cityblockBox = tk.StringVar()

			#create widgets
			self.EnsemSecondWindow = ttk.Label(self,text="EnsembleClustring").grid(column=1,row=0,columnspan=3,sticky=(N))
			self.DistMetInterest = ttk.Label(self,text="Which distance metrics would you like to use?").grid(column=1,row=1,columnspan=3,sticky=(N))
			self.euclideanCheck = ttk.Checkbutton(self,text="Euclidean",variable=euclideanBox,onvalue='euclidean',offvalue='Not').grid(column=1,row=2,sticky=(N,W))
			self.seuclideanCheck = ttk.Checkbutton(self,text="Standardized Euclidean",variable=seuclideanBox,onvalue='seuclidean',offvalue='Not').grid(column=2,row=2,sticky=(N,W))
			self.sqeuclideanCheck = ttk.Checkbutton(self,text="Square Root-Euclidean",variable=sqeuclideanBox,onvalue='sqeuclidean',offvalue='Not').grid(column=3,row=2,sticky=(N,W))
			self.cosineCheck = ttk.Checkbutton(self,text="Cosine",variable=cosineBox,onvalue='cosine',offvalue='Not').grid(column=1,row=3,sticky=(N,W))
			self.chebyshevCheck = ttk.Checkbutton(self,text="Chebyshev",variable=chebyshevBox,onvalue='chebyshev',offvalue='Not').grid(column=2,row=3,sticky=(N,W))
			self.correlationCheck = ttk.Checkbutton(self,text="Correlation",variable=correlationBox,onvalue='correlation',offvalue='Not').grid(column=3,row=3,sticky=(N,W))
			self.canberraCheck = ttk.Checkbutton(self,text="Canberra",variable=canberraBox,onvalue='canberra',offvalue='Not').grid(column=1,row=4,sticky=(N,W))
			self.braycurtisCheck = ttk.Checkbutton(self,text="Bray-Curtis",variable=braycurtisBox,onvalue='braycurtis',offvalue='Not').grid(column=2,row=4,sticky=(N,W))
			self.minkowskiCheck = ttk.Checkbutton(self,text="Minkowski",variable=minkowskiBox,onvalue='minkowski',offvalue='Not').grid(column=3,row=4,sticky=(N,W))
			self.cityblockCheck = ttk.Checkbutton(self,text="City Block",variable=cityblockBox,onvalue='cityblock',offvalue='Not').grid(column=1,row=5,sticky=(N,W))
			self.SecondNextButton = ttk.Button(self, text="Next->", command = secondNext).grid(column=2,row=6,sticky=(N))
			self.home1 = ttk.Button(self, text="Return to Home",command=self.home).grid(column=2,row=7,sticky=(N))

		def secondNext(*args):
			global distanceMetList
			euclideanBoxOut = euclideanBox.get()
			seuclideanBoxOut = seuclideanBox.get()
			sqeuclideanBoxOut = sqeuclideanBox.get()
			cosineBoxOut = cosineBox.get()
			chebyshevBoxOut = chebyshevBox.get()
			correlationBoxOut= correlationBox.get()
			canberraBoxOut = canberraBox.get()
			braycurtisBoxOut = braycurtisBox.get()
			minkowskiBoxOut = minkowskiBox.get()
			cityblockBoxOut = cityblockBox.get()


			#input distances measures check values into a list and determine if they have been selected. 
			distMetL = [euclideanBoxOut,seuclideanBoxOut,sqeuclideanBoxOut,cosineBoxOut,chebyshevBoxOut,correlationBoxOut,canberraBoxOut,braycurtisBoxOut,minkowskiBoxOut,cityblockBoxOut]
			distanceMetList = []
			for i in range(len(distMetL)):
				#distance metric list for usage
				if len(distMetL[i]) > 0:
					distanceMetList.append(distMetL[i])
			
			#make linkParams global for easy access when sending to ensemble clustering.
			global linkParams
			linkParams =[]
			if 'ward' in ensemLinkList:
				#put the ward euclidean measure in and remove ward from the list
				ensemLinkList.remove('ward')
				linkParams.append(['ward','euclidean'])
				for i in range(len(ensemLinkList)):
					for k in range(len(distanceMetList)):
						linkParams.append([ensemLinkList[i], distanceMetList[k]])

			else:
				#put all the linkage parameters into list of lists for ensemble clustering.
				for i in range(len(ensemLinkList)):
					for k in range(len(distanceMetList)):
						linkParams.append([ensemLinkList[i], distanceMetList[k]])

			### Selection Error handling------------------------------------------------------------------------------------------------------------------------------------------------
			if len(linkParams) == 1:
				messagebox.showerror(title="Selection Error", message="You cannot run an ensemble, you have only selected one linkage-distance pair. Please select at least one other distance metric.")
				return
			###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			objects = self.grid_slaves()
			for i in objects:
				i.grid_forget()
			
			#create widgets for the clustergram function input. 
			self.EnsembleLabel = ttk.Label(self, text="Ensemble Clustering",font=("TkHeadingFont",24)).grid(column=0,row=0,sticky=(N),columnspan=3)
			self.JuneLab = ttk.Label(self, text="# of Clusters",font=("TkHeadingFont",18)).grid(column=0,row=4,sticky=(N))
			self.MetabNumPClustL = ttk.Label(self, text="Min. # of features?", font=("TkHeadingFont",18)).grid(column=1,row=4,sticky=(N))
			self.colMap = ttk.Label(self, text="ColorMap", font=("TkHeadingFont",18)).grid(column=0,row=1,sticky=(N))
			self.trans = ttk.Label(self,text="Transform",font=("TkHeadingFont",18)).grid(column=1,row=1,sticky=(N))
			self.scale = ttk.Label(self,text="Scale", font=("TkHeadingFont",18)).grid(column=2,row=1,sticky=(N))
			self.mstF1 = ttk.Button(self,text="Run MST",command=self.mstF).grid(column=0,row=8,sticky=(N))
			self.ButtonColormaps = ttk.Button(self,text="ColorMap Options",command=cMapOpt).grid(column=0,row=3,sticky=(N))
			self.ButtonTransforms = ttk.Button(self,text="Transform Opttions",command=transOpts).grid(column=1,row=3,sticky=(N))
			self.ButtonScale = ttk.Button(self,text="Scale Options", command=scaleOpts).grid(column=2,row=3,sticky=(N))
			self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=2,row=8, sticky=(N))
			
			global colorList
			colorList = config.colorList
			
			#create list box of the ensemble optimal clusters
			global optClustBox
			global minNumMetabClust
			global colorMapEnsem
			global transformEnsem
			global scaleEnsem
			optClustBox = Listbox(self,height=5,width=35)
			minNumMetabClust = Listbox(self, height=5,width=35)
			colorMapEnsem = Listbox(self, height=5,width=35)
			transformEnsem = Listbox(self, height=5,width=35)
			scaleEnsem = Listbox(self, height=5,width=35)

			global transformList
			global scaleList
			transformList = config.transformList
			scaleList = config.scaleList


			#Create the lists of available options for selection 
			optClusters = tuple(range(1,101))
			minMetab = tuple(range(0,101))
			distNames = StringVar(value=optClusters)
			minNum = StringVar(value=minMetab)
			colMap = StringVar(value=colorList)
			transType = StringVar(value=transformList)
			scaleType = StringVar(value=scaleList)

			for i in range(len(colorList)):
				colorMapEnsem.insert(i,colorList[i])
			
			#create binding event and 
			optClustBox.bind('<Double-1>',ensembleOptClust)
			scaleEnsem.bind('<Double-1>',ensembleScale)
			colorMapEnsem.bind('<Double-1>',ensembleCMap)
			transformEnsem.bind('<Double-1>',ensembleTransform)
			self.optClustBox = optClustBox
			self.optClustBox.grid(column=0,row=5,columnspan=1)
			self.minNumMetabClust = minNumMetabClust
			self.minNumMetabClust.grid(column=1,row=5,columnspan=1)
			self.colorMapEnsem = colorMapEnsem
			self.colorMapEnsem.grid(column=0,row=2,columnspan=1)
			self.transformEnsem = transformEnsem
			self.transformEnsem.grid(column=1,row=2,columnspan=1)
			self.scaleEnsem = scaleEnsem
			self.scaleEnsem.grid(column=2,row=2,columnspan=1)

			
		def ensembleOptClust(*args):
			#grab the current selection of the list
			global optClust
			optClust = optClustBox.curselection()
			optClust = int(optClust[0])+1
			minMetab = tuple(range(0,101))

			lenList = len(minNumMetabClust.get(0,tk.END))
			if lenList > 0:
				minNumMetabClust.delete(0,lenList-1)

			for i in range(len(minMetab)):
				minNumMetabClust.insert(i,minMetab[i])
			self.submit = ttk.Button(self,text="Submit",command=ensembleGo).grid(column=1,row=8,sticky=(N))
		
		def ensembleScale(*args):
			#grab the selection and minimum number of clusters
			global scaleSel
			scaleSel = scaleEnsem.curselection()
			scaleSel = scaleList[scaleSel[0]]
			optClusters = tuple(range(1,101))
			lenList = len(optClustBox.get(0,tk.END))
			if lenList > 0:
				optClustBox.delete(0,lenList-1)

			for i in range(len(optClusters)):
				optClustBox.insert(i,optClusters[i])

		def ensembleCMap(*args):
			#grab the selection and colormap wanted
			global cMapSel
			cMapSel = colorMapEnsem.curselection()
			cMapSel = colorList[cMapSel[0]]
			lenList = len(transformEnsem.get(0,tk.END))
			if lenList > 0:
				transformEnsem.delete(0,lenList-1)

			for i in range(len(transformList)):
				transformEnsem.insert(i,transformList[i])

		def ensembleTransform(*args):
			#grab the data transform of interest
			global transSel
			transSel = transformEnsem.curselection()
			transSel = transformList[transSel[0]]
			lenList = len(scaleEnsem.get(0,tk.END))
			if lenList > 0:
				scaleEnsem.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleEnsem.insert(i,scaleList[i])

		def ensembleGo(*args):
			minFeature = minNumMetabClust.curselection()
			minFeature = int(minFeature[0])


			#send to ensemble clustering algorithm
			GU.ensembleClustering(optNum= optClust, minMetabs= minFeature, colorMap = cMapSel,linkParams=linkParams,transform=transSel,scale=scaleSel)

		def cMapOpt(*args):
			#send to webpage of colormap options in python
			webbrowser.open('https://matplotlib.org/stable/tutorials/colors/colormaps.html')

		def transOpts(*args):
			webbrowser.open('https://github.com/hisl6802/Transformation-and-Scaling/wiki/Transformations')

		def scaleOpts(*args):
			webbrowser.open('https://github.com/hisl6802/Transformation-and-Scaling/wiki/Scaling')

		def standardEnsemble(*args):
			objects = self.grid_slaves()
			for i in objects:
				i.grid_forget()
			
			#set linkage params as a global for ease of access when selecting standard ensemble.
			global linkParams
			linkParams = [['ward','euclidean'],['single','euclidean'],['single','sqeuclidean'],['single','seuclidean'],['single','chebyshev'],['complete','euclidean'],['complete','sqeuclidean'],['complete','seuclidean'],['complete','chebyshev'],['average','euclidean'],['average','sqeuclidean'],['average','seuclidean'],['average','chebyshev']]
			
			#create widgets for the clustergram function input. 
			self.EnsembleLabel = ttk.Label(self, text="Ensemble Clustering",font=("TkHeadingFont",24)).grid(column=0,row=0,sticky=(N),columnspan=3)
			self.JuneLab = ttk.Label(self, text="# of Clusters",font=("TkHeadingFont",18)).grid(column=0,row=4,sticky=(N))
			self.MetabNumPClustL = ttk.Label(self, text="Min. # of features?", font=("TkHeadingFont",18)).grid(column=1,row=4,sticky=(N))
			self.colMap = ttk.Label(self, text="ColorMap", font=("TkHeadingFont",18)).grid(column=0,row=1,sticky=(N))
			self.trans = ttk.Label(self,text="Transform",font=("TkHeadingFont",18)).grid(column=1,row=1,sticky=(N))
			self.scale = ttk.Label(self,text="Scale", font=("TkHeadingFont",18)).grid(column=2,row=1,sticky=(N))
			self.ButtonColormaps = ttk.Button(self,text="ColorMap Options",command=cMapOpt).grid(column=0,row=3,sticky=(N))
			self.ButtonTransforms = ttk.Button(self,text="Transform Options",command=transOpts).grid(column=1,row=3,sticky=(N))
			self.ButtonScale = ttk.Button(self,text="Scale Options", command=scaleOpts).grid(column=2,row=3,sticky=(N))
			self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=2,row=8, sticky=(N))
			
			global colorList
			colorList = config.colorList 
			
			#create list box of the ensemble optimal clusters
			global optClustBox
			global minNumMetabClust
			global colorMapEnsem
			global transformEnsem
			global scaleEnsem
			optClustBox = Listbox(self,height=5,width=35)
			minNumMetabClust = Listbox(self, height=5,width=35)
			colorMapEnsem = Listbox(self, height=5,width=35)
			transformEnsem = Listbox(self, height=5,width=35)
			scaleEnsem = Listbox(self, height=5,width=35)

			global transformList
			global scaleList
			transformList = config.transformList
			scaleList = config.scaleList


			#Create the lists of available options for selection 
			optClusters = tuple(range(1,101))
			minMetab = tuple(range(0,101))
			distNames = StringVar(value=optClusters)
			minNum = StringVar(value=minMetab)
			colMap = StringVar(value=colorList)
			transType = StringVar(value=transformList)
			scaleType = StringVar(value=scaleList)

			for i in range(len(colorList)):
				colorMapEnsem.insert(i,colorList[i])
			
			#create binding event and 
			optClustBox.bind('<Double-1>',ensembleOptClust)
			scaleEnsem.bind('<Double-1>',ensembleScale)
			colorMapEnsem.bind('<Double-1>',ensembleCMap)
			transformEnsem.bind('<Double-1>',ensembleTransform)
			self.optClustBox = optClustBox
			self.optClustBox.grid(column=0,row=5,columnspan=1)
			self.minNumMetabClust = minNumMetabClust
			self.minNumMetabClust.grid(column=1,row=5,columnspan=1)
			self.colorMapEnsem = colorMapEnsem
			self.colorMapEnsem.grid(column=0,row=2,columnspan=1)
			self.transformEnsem = transformEnsem
			self.transformEnsem.grid(column=1,row=2,columnspan=1)
			self.scaleEnsem = scaleEnsem
			self.scaleEnsem.grid(column=2,row=2,columnspan=1)

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create first window
		self.EnsembleLabel = ttk.Label(self,text="Ensemble Clustering",font=("TkHeadingFont",18)).grid(column=0,row=0,sticky=(N))
		self.Win1Q = ttk.Label(self,text="Which linkage functions would you like to include?",font=("TkHeadingFont",16)).grid(column=0,row=1,sticky=(N))
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=7, sticky=(N),columnspan=2)
		self.FirstNext = ttk.Button(self,text="Next->",command=firstNext).grid(column=0,row=6,sticky=(N))
		self.standardEnsemble = ttk.Button(self, text="Pre-set Ensemble",command=standardEnsemble).grid(column=0,row=8)
		singleBox = tk.StringVar()
		completeBox = tk.StringVar()
		averageBox = tk.StringVar()
		wardBox = tk.StringVar()
		self.SingleCheck = ttk.Checkbutton(self,text="Single",variable=singleBox,onvalue='single',offvalue='Not').grid(column=0,row=2,sticky=(N))
		self.CompleteCheck = ttk.Checkbutton(self,text="Complete",variable=completeBox,onvalue='complete',offvalue='Not').grid(column=0,row=3,sticky=(N))
		self.AverageCheck = ttk.Checkbutton(self,text="Average",variable=averageBox,onvalue='average',offvalue='Not').grid(column=0,row=4,sticky=(N))
		self.WardCheck = ttk.Checkbutton(self,text="Ward",variable=wardBox,onvalue='ward',offvalue='Not').grid(column=0,row=5,sticky=(N))


	def mst(self):
		def valType(*args):
			selection = valTypeBox.curselection()
			index = selection[0]
			selection =valList[selection[0]]
			#send the user to the minimum spanning tree function.
			GU.MST(self, transform = transform, scale=scale,func=selection)

		def transType(*args):
			global transform
			transform = transBox.curselection()
			transform = transformList[transform[0]]
		
			lenList = len(scaleBox.get(0,tk.END))
			if lenList > 0:
				scaleBox.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleBox.insert(i,scaleList[i])

		def scaleType(*args):
			global scale
			scale = scaleBox.curselection()
			scale = scaleList[scale[0]]

			lenList = len(valTypeBox.get(0,tk.END))
			if lenList > 0:
				valTypeBox.delete(0,lenList-1)

			for i in range(len(valList)):
				valTypeBox.insert(i,valList[i])


		#delete previous objects on the GUI
		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Validation Measure",font=("TkHeadingFont",24)).grid(column=1,row=0,sticky=(N))
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=1,row=3, sticky=(N))
		self.transLab =ttk.Label(self,text="Transform",font=("TkHeadingFont",16)).grid(column=0,row=1,sticky=(N))
		self.scaleLab =ttk.Label(self,text="Scale",font=("TkHeadingFont",16)).grid(column=1,row=1,sticky=(N))
		self.valLab = ttk.Label(self,text="Validation",font=("TkHeadingFont",16)).grid(column=2,row=1,sticky=(N))

		#validation index list (MST-based, DBI, Dunn)
		transformList = config.transformList
		scaleList = config.scaleList
		valList = ('k-means based','DBI','Dunn','PBM','Silhouette')
		valTypeBox = Listbox(self,height=5,width=30)
		transBox = Listbox(self,height=5,width=30)
		scaleBox = Listbox(self,height=5,width=30)

		#defining string variables for each listbox
		transName = StringVar(value=transformList)
		scaleName = StringVar(value=scaleList)
		valName = StringVar(value=valTypeBox)

		for i in range(len(transformList)):
			transBox.insert(i, transformList[i])

		#create a binding event
		transBox.bind('<Double-1>',transType)
		scaleBox.bind('<Double-1>',scaleType)
		valTypeBox.bind('<Double-1>',valType)
		self.transBox = transBox
		self.transBox.grid(column=0,row=2,sticky=(N),padx=5,pady=5)
		self.scaleBox =scaleBox
		self.scaleBox.grid(column=1,row=2,sticky=(N),pady=5)
		self.valTypeBox = valTypeBox
		self.valTypeBox.grid(column=2,row=2,sticky=(N),padx=5,pady=5)

	def genSelClustFig(self):

		def selectedCMap(*args):
			clMap = cmap.curselection()
			clMap = colorList[clMap[0]]
			GB.createHeatmapFig(clMap=clMap)

		objects =self.grid_slaves()
		for i in objects:
			i.grid_forget()

		self.Label = ttk.Label(self,text='Select ColorMap', font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N))
		cmap1= tk.StringVar()
		cmap=Listbox(self,width=25)
		self.home2 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=1,row=2,sticky=(N))
		colorList = config.colorList

		cmap.bind('<Double-1>',selectedCMap)

		for i in range(len(colorList)):
			cmap.insert(i,colorList[i])

		self.cmap = cmap
		self.cmap.grid(column=1,row=1,sticky=(N))


if __name__ == '__main__':
	#launch application
	root = tk.Tk()

	#ask for the directory the user would like to save the files to. 
	messagebox.showinfo(title="INFO",message="Please select directory where you would like the log and output files to be saved!")
	#get directory where user would like to save the log and output files.
	directory = filedialog.askdirectory()
	os.chdir(directory)

	app = JuneLabClusteringGUI(master=None)
	app.master.minsize(200,200)
	height = app.master.winfo_screenheight() - 400
	width = (app.master.winfo_screenwidth())*0.9
	width = int(width)
	app.master.maxsize(900,900)
	width = str(width)
	height = str(height)
	screenSize = width +'x' + height + '+0+100'
	app.master.geometry("")
	app.master.resizable(0,0)

	app.master.title("Home")
	app.mainloop()
