# Directory to the mesh
set dir /home/mroe734/Documents/phd/mesh/

# Read elements and nodes of the volumetric mesh
gfx read node AWA015_PTA_1_Rec_Trans_volumetric_mesh_annotated.exnode 
gfx read elem AWA015_PTA_1_Rec_Trans_volumetric_mesh_annotated.exelem
gfx def faces egroup uterus

# Create scene and spectrum
gfx create win
gfx create spectrum thickness

# Create uterus mesh and colour it
gfx modify g_element /uterus surfaces data thickness spectrum thickness
gfx modify spectrum thickness autorange

# Create colour bar
gfx create colour_bar spectrum thickness
gfx modify g_element /uterus point glyph colour_bar spectrum thickness LOCAL NORMALISED_WINDOW_FIT_LEFT

gfx modify window 1 image view_all

