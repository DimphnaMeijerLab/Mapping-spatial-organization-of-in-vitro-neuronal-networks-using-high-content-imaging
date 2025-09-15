//-----------------------------------------------------------------------
// Show Tile Numbers of Stitched Images
// Macro for Showing the tiles of stitched images (from high content microscope)
// Created 20-06-2022 by Martijn Joosten.
//
// 
//
//
//---------------------------START FUNCTIONS-----------------------------


function initialize_position(width, height, clockwise){

	//-----------------------------------------------------
	// Finds the center of a spiral grid, i.e. the
	// (x,y) location in the spiral grid with value 0.

	// Inputs
	// width & height: the dimensions of the spiral grid.

	// Output
	// init_position: array of length 2.
	// Index 0 is the x position, index 1 is the y position.
	//-----------------------------------------------------
	add_y = 0;
	if(clockwise == false){
		add_y = 1;
	}

	init_position = newArray(2);
	if (width%2 == 0){
		x = width / 2 - 1;
	}else{
		x = floor(width / 2);
	}
	if (height%2 == 0){
		y = height / 2 - 1 + add_y;
	}else{
		y = floor(height / 2);
	}
	init_position[0] = x;
	init_position[1] = y;

	return init_position;

}

function compare_arrays(arr1, arr2){

	//----------------------------------------------------------
	// Function compares array arr1 and array arr2 and
	// returns true if they are equal, and false if they are not.
	//----------------------------------------------------------

	same = true;
	l = arr1.length;
	for (i = 0; i < l; i++){
		if (arr1[i] != arr2[i]){
			same = false;
		}
	}
	return same;
}

function turn_right(current_direction, clockwise){

	//-----------------------------------------------------
	// This function is called by make_spiral_grid().

	// It takes as input the current direction (north, south, west and east),
	// and outputs the direction after turning right:
	// if clockwise:
	// north->east, south->west, east->south, west->north.
	// else (counterclockwise):
	// north->west, south->east, east->north, west->south.

	// All directions are arrays of length 2:
	// index 0 is dx, index 1 is dy.
	//-----------------------------------------------------

	NORTH = newArray(0,-1);
	S = newArray(0,1);
	W = newArray(-1,0);
	E = newArray(1,0);

	if (clockwise){
		if (compare_arrays(current_direction, NORTH)){
			new_direction = E;
		}else if (compare_arrays(current_direction, S)){
			new_direction = W;
		}else if (compare_arrays(current_direction, E)){
			new_direction = S;
		}else if (compare_arrays(current_direction, W)){
			new_direction = NORTH;
		}

	}else{
		if (compare_arrays(current_direction, NORTH)){
			new_direction = W;
		}else if (compare_arrays(current_direction, S)){
			new_direction = E;
		}else if (compare_arrays(current_direction, E)){
			new_direction = NORTH;
		}else if (compare_arrays(current_direction, W)){
			new_direction = S;
		}
	}

	return new_direction;

}

function make_spiral_grid(width,height,clockwise){
	//-----------------------------------------------------
	// This function makes a spiral grid array.
	// Note that 2D arrays are not supported by Fiji macro language.
	// The output will be a 1D array with length width*height.

	// If you index it with matrix[x + y * width],
	// the matrix is effectively 2D.
	//-----------------------------------------------------

	// Check if dimenstions are at least 1:
	if (width < 1 || height < 1){
		return;
	}

	// Directions to walk (dx, dy) are arrays:
	NORTH = newArray(0,-1);
	S = newArray(0,1);
	W = newArray(-1,0);
	E = newArray(1,0);

	// Initial position:
	init_position = initialize_position(width, height, clockwise);
	x = init_position[0];
	y = init_position[1];

	// We want to start walking to the west.
	// This means our initial direction is north (clockwise) or south (counterclockwise):
	// we then turn right immediately, and end up going west.
	if(clockwise){
		direction = NORTH;
	}else{
		direction = S;
	}
	dx = direction[0];
	dy = direction[1];

	// Initialize matrix:
	matrix = newArray(width*height);
	Array.fill(matrix, NaN);
	count = 0;

	while (true){
		// Fill matrix at current position with value count:
		matrix[x + y*width] = count;

		count = count + 1;

		// Try to turn right:
		new_direction = turn_right(direction, clockwise);
		new_dx = new_direction[0];
		new_dy = new_direction[1];
		new_x = x + new_dx;
		new_y = y + new_dy;

		// Turn right if you are not yet at the boundaries,
		// and if the position at your right-hand side is not yet visited:
		if ((0<=new_x) & (new_x<width) & (0<=new_y) & (new_y<height) & (isNaN(matrix[new_x + new_y*width]))){
			x = new_x;
			y = new_y;
			dx = new_dx;
			dy = new_dy;
			direction = new_direction;
		}
		// If not, go straight:
		else{
			x = x + dx;
			y = y + dy;
			// If you are at the boundary, stop the process and return the matrix:
			if (!(((0 <= x) & (x < width)) & (0 <= y) & (y < height))){
				return matrix;
			}
		}
	}
}

function make_column_grid(width, height){
	//-----------------------------------------------------
	// This function makes a column grid array.
	// Note that 2D arrays are not supported by Fiji macro language.
	// The output will be a 1D array with length width*height.

	// If you index it with matrix[x + y * width],
	// the matrix is effectively 2D.
	//-----------------------------------------------------
	matrix = newArray(width*height);
	Array.fill(matrix, NaN);
	count = 0;
	x = 0;
	y = 0;
	for (i = 0; i < width*height; i++) {
		matrix[x + y*width] = count;
		count = count + 1;
		y = y + 1;
		if (y == height){
			y = 0;
			x = x + 1;
		}
	}
	return matrix;
}



choice = newArray("Top", "Left", "Center", "Right", "Bottom");

Dialog.createNonBlocking("Alignment");
Dialog.addNumber("size1", 15);
Dialog.addNumber("size2", 15);
Dialog.addNumber("w", 11);
Dialog.addNumber("h", 11);
Dialog.addChoice("Select Alignment", choice, "Yes");
Dialog.show();

size1 = Dialog.getNumber();
size2 = Dialog.getNumber();
w = Dialog.getNumber();
h = Dialog.getNumber();

Alignment = Dialog.getChoice();

column_grid = make_column_grid(size1,size2);
spiral_grid = make_spiral_grid(size1,size2,true);



diff1 = (size1 - h)/2;
diff2 = (size2 - w)/2;

W = getWidth();
H = getHeight();
width = round(W / w);
height = round(H / h);
setJustification("center");
setFont("SansSerif", 200);
setColor("white");

x = 0;
y = 0;

tile_nr_column=newArray();
tile_nr_spiral=newArray();


//Get tile numbers for the column grid
for (i=0; i<size1; i++) {
	for (j=0; j<size2; j++) {
		if (Alignment == "Top") {
			if (i>=diff1 && i<size1-diff1) {
				if (j>=0 && j<size2-diff2*2) {
					tile_nr_column = Array.concat(tile_nr_column, column_grid[i + (j+1)*size1]);
					tile_nr_spiral = Array.concat(tile_nr_spiral, spiral_grid[i + (j+1)*size1]);
				}
			}
		}
		if (Alignment == "Left") {
			if (i>=0 && i<size1-diff2*2) {
				if (j>=diff2 && j<size2-diff2) {
					tile_nr_column = Array.concat(tile_nr_column, column_grid[i+1 + (j)*size1]);
					tile_nr_spiral = Array.concat(tile_nr_spiral, spiral_grid[i+1 + (j)*size1]);
				}
			}
		}
		if (Alignment == "Center") {
			if (i>=diff1 && i<size1-diff1) {
				if (j>=diff2 && j<size2-diff2) {
					tile_nr_column = Array.concat(tile_nr_column, column_grid[i + (j)*size1]);
					tile_nr_spiral = Array.concat(tile_nr_spiral, spiral_grid[i + (j)*size1]);
				}
			}
		}
		if (Alignment == "Right") {
			if (i>=diff1*2 && i<size1) {
				if (j>=diff2 && j<size2-diff2) {
					tile_nr_column = Array.concat(tile_nr_column, column_grid[i-1 + (j)*size1]);
					tile_nr_spiral = Array.concat(tile_nr_spiral, spiral_grid[i-1 + (j)*size1]);
				}
			}
		}
		if (Alignment == "Bottom") {
			if (i>=diff1 && i<size1-diff1) {
				if (j>=diff2*2 && j<size2) {
					tile_nr_column = Array.concat(tile_nr_column, column_grid[i + (j-1)*size1]);
					tile_nr_spiral = Array.concat(tile_nr_spiral, spiral_grid[i + (j-1)*size1]);
				}
			}
		}
	}
}


y=0;
x=0;
//Draw the tile numbers for column and spiral grid over the image
for (i = 0; i < w*h; i++) {
	drawString(tile_nr_spiral[i], x*width + width/2, y*height + height/2);
	drawString(tile_nr_column[i], x*width + width/2, y*height + height/2+200);
	// Update position
	y = y + 1;
	if (y == h){
		y = 0;
		x = x + 1;
	}
}
