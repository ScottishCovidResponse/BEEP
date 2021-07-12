// Lists all the variables used in the graphical interface

var	width, height, heightpage, widthold, heightold;       // Define the size of the window

var initdone = 0;                                         // This is set to one when the initialisation is done

var rateradio = "rate";                                   // Radio button for if rate is displayed
var lowerhighlight = "off";                               // Dertermines if a high is added

var visjson;                                              // This object stores information about graphs

var plottype;                                             // Sets the plot type

var zoomfac, zoomx, zoomy;                                // Sets up the zoom for maps

var menuslidest, menuslidesize;                           // Information for menu slider

var tree=[];                                              // Stores the tree structure for the menu

var col_map=[];                                           // Creates a colour map for displaying map 

var selparam_table;                                       // The table for a selected parameter
var selparam_name;     

var seleq;                                                // The list of parameters for a selected equation
var seleq_name;                                           // The name of the selected equation

var paramsel;                                             // The name of a selected parameter

var check;                                                // Stores the status of the check button                                        

var param_link=[];                                        // Provides link to the model 

var playtime = 0;                                         // When performing animations this is the plot time
var playing = 0, playtimemax, playstartt, playtimeinit;   // Used for playing animations

var slidey = 0, slidefrac, slidey1, slidey2;              // Used for the slider on the right menu

var tabx_slide = 0, tabx_slidefrac, tabx1, tabx2;         // Used for the x slider on tables
var taby_slide = 0, taby_slidefrac, taby1, taby2;         // Used for the y slider on tables

var clickline = 100;                                      // The number of horizontal lines used to make clickable map
var areaover = -1;                                        // Keeps track of being over an area

var sourceon = 0;                                         // Determines if source is shown

var maincv, graphcv, cv, resultcv;                        // Used for drawing on canvases in HTML5

var pagesub=[], pagesubsub=[], pagesubsubsub=[];          // Stores the page/subpage being viewed
var page=0, ps=0, pss=0, psss=0;                          // Stores the current page position

var over = -1;                                            // The button the mouse is over

var model;                                                // Stores the compartmental model

var foi;                                                  // Stores the equation for the force of infection
 
var canover = -1;                                         // The canvas button the mouse is over

var nticksize, ticksize=[];                               // Used to determine ticks along axes

var ntickx, tickx=[], nticky, ticky=[];                   // Stores information about ticks on the axes

var axxmin, axxmax, axymin, axymax;                       // The range along the x and y axes

var nbut, buttext=[], butx=[], buty=[];                   // Stores information about buttons
var butdx=[], butdy=[], butac=[], buttype=[];
var butover=[], butval=[], butval2=[];

var ncanbut, canbuttext=[], canbutx=[], canbuty=[];       // Stores information about canvas buttons
var canbutdx=[], canbutdy=[], canbutac=[], canbuttype=[];
var canbutover=[], canbutval=[], canbutval2=[];

var polypoint=[];                                         // Used for plotting polygons
for(i = 0; i < 10; i++) polypoint[i] = []; 

var nlines, lines=[], linesy=[], lineheight = 23;         // Used when drawing paragraphs

var timedown, timeup;                                     // The timing the mouse was pushed down and went up

var arrow = 0;                                            // Determines the pointer (arrow or hand)

var mx, my;                                               // The mouse position

var mxst, myst, drag = 0;                                 // Used when dragging objects

var alertdone = 0;                                        // This is set to 1 if an alert has been called
