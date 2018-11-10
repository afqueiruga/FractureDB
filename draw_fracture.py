from PIL import Image, ImageDraw

def draw_fracture(x1,y1,x2,y2, size=(31,31)):
    im = Image.new('1', size, color=1)
    xy_to_px = lambda x,y: ( (x+1)/2.0*size[0],
                             (y+1)/2.0*size[1] )
    draw = ImageDraw.Draw(im)
    draw.line( xy_to_px(-1,0) + xy_to_px(x1,y1), fill=0)
    draw.line( xy_to_px(x1,y1) + xy_to_px(x2,y2), fill=0)
    del draw
    #im.save(sys.stdout, "PNG")
    return im
