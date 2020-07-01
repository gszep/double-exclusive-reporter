data_path = './20191121.tif'
name,ext = splitext(basename(data_path))

grid_height,grid_width = 5,4
area_threshold = 1e5

hyperstack = io.imread(data_path)
ntime,height,width,channels = hyperstack.shape

assert channels == 3
cfp,yfp,rfp = 0,1,2

gridmask = hyperstack[-1,:,:,rfp] > mean(hyperstack[-1,:,:,rfp])
gridmask = binary_closing(binary_opening(gridmask,structure=ones((10,10))),structure=ones((20,20)))

region_boxes =  array([ region.bbox for region in regionprops(label(clear_border(gridmask))) if region.area < area_threshold ])
region_centre = array([ [(box[0]+box[2])/2,(box[1]+box[3])/2] for box in region_boxes  ])

region_index =  array([ [int(grid_height*y/height)  ,int(grid_width*x/width)   ] for y,x in region_centre ])
regions = list(zip(region_index,region_boxes))

fig = figure(figsize=(12,12))
imshow(hyperstack[-1,:,:,yfp],cmap='yellow')
imshow(hyperstack[-1,:,:,cfp],cmap='cyan')

for (i,j),(miny, minx, maxy, maxx) in regions:
    text(minx,miny,  '{},{}'.format(i,j))

    hyperstack_region = hyperstack[:,miny:maxy,minx:maxx,:]
    io.imsave(join(name,'{}-{}.tif'.format(i,j)),hyperstack_region)

    rectangle = Rectangle((minx,miny), maxx-minx, maxy-miny)
    fig.axes[0].add_collection(PatchCollection([rectangle],facecolor='',edgecolor='red'))
