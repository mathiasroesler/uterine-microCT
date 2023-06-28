set dir /home/mroe734/Documents/phd/microCT/data/AWA015_PTA_1_Rec_Trans/downsampled/ST/binary

# Read the data
gfx read node Streamlines_L4_FB
gfx read elem Streamlines_L4_FB

# Create the scene
gfx create win

# Create the spectrum
gfx create spectrum angle linear reverse range 0 90

# Create fibre lines and colour them
gfx modify g_element /streamlines lines data angle spectrum angle line_width 2

# Create colour bar
gfx create colour_bar spectrum angle number_format %.1e
gfx modify g_element /streamlines point glyph colour_bar spectrum angle LOCAL NORMALISED_WINDOW_FIT_LEFT

# Display everything
gfx modify window 1 image view_all
