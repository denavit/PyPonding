from structures import ExampleRoof

roof = ExampleRoof()
roof.use_CBDI = False

#roof.plot_load_cells = True

roof.BuildModel()
#roof.PlotModel()
print(roof.ColumnReaction('B2'))
print(roof.ColumnReaction('C2'))
