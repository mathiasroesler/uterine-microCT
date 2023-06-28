set dir /home/mroe734/Documents/phd/microCT/data/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/binary

# Read the data
gfx read node Streamlines_L4_FB
gfx read elem Streamlines_L4_FB

# Create the scene
gfx create win
gfx define font large "21 decorative normal bold"

# Create the spectrum
gfx create spectrum angle linear reverse range 0 90

# Create fibre lines and colour them
gfx modify g_element /streamlines lines data angle spectrum angle line_width 2

# Create colour bar
gfx create colour_bar spectrum angle number_format %.1e font large centre -1.8 0 0.5
gfx modify g_element /streamlines point glyph colour_bar spectrum angle LOCAL NORMALISED_WINDOW_FIT_LEFT scale_factors 0.5

# Display everything
gfx modify window 1 image view_all

# Rotate the image
gfx modify window 1 image rotate 1 0 0 -90
gfx modify window 1 image rotate 0 1 0 45
#gfx modify window 1 image rotate 0 0 1 -180

gfx modify window 1 layout two_free
gfx modify window 1 image rotate 0 0 1 90
gfx modify window 1 image rotate 1 0 0 180

