/* 
 * Image3D2ViewConfig.java
 *
 * Class for holding attributes that determine various visual aspects of
 *      the geometric 3D (ie Image3D2) AceTree window.
 *
 * 2.14.2014
 */

package org.rhwlab.image;

import java.awt.Color;

public class Image3DViewConfig {

    ////////////////////////////////////////////////////////////
    // The following attributes are for various aspects of the
    // Color Controls panel
    ////////////////////////////////////////////////////////////

    //////////////////// DEFAULTS //////////////////// 

    // color choices that the user can assign to a given lineage
    public final static String [] LINEAGE_COLORS = {
        "red"
        ,"blue"
        ,"green"
        ,"yellow"
        ,"cyan"
        ,"magenta"
        ,"pink"
        ,"gray"
        ,"light gray"
        ,"dark gray"
        ,"white"
        ,"transparent"
        ,"omit"
    };

    // given a String for a color, sees if it is in LINEAGE_COLORS
    // if it is, returns its index
    // else returns 0
    public static int getColorIndex(String color) {
        for(int i = 0; i < LINEAGE_COLORS.length; i++)
            if(color.equals(LINEAGE_COLORS[i]))
                return i;
        return 0;
    }

    // used for 'other' option but not really sure what good it is
    // seems to have been joined with LINEAGE_COLORS so will consider removing
    public final static String [] TRANSPROPS = {
        "omit"
        ,"transparent"
        ,"white"
    };

    // color options for background color of main image display
    public final static String [] GRAYDEPTH = {
        "white"
        ,"light gray"
        ,"dark gray"
    };

    // represents the sublineage and associated color 
    // displayed in Lineage Color Controls
    // numbers correspond to index of color in LINEAGE_COLORS
    //          note: add button to add/remove display prop
    public final SublineageDisplayProperty [] DISP_PROPS = {
         new SublineageDisplayProperty("ABa", 0)
	    ,new SublineageDisplayProperty("ABp", 1)
	    ,new SublineageDisplayProperty("C", 5)
	    ,new SublineageDisplayProperty("D", 6)
	    ,new SublineageDisplayProperty("E", 2)
	    ,new SublineageDisplayProperty("MS", 4)
	    ,new SublineageDisplayProperty("P", 3)
	    ,new SublineageDisplayProperty("polar", 7)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
        ,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("other", 2)
	    ,new SublineageDisplayProperty("background", 1)
    };
    
    public static final SublineageDisplayProperty [] DEFAULT_DISP_PROPS = {
    	 new SublineageDisplayProperty("ABa", 0)
	    ,new SublineageDisplayProperty("ABp", 1)
	    ,new SublineageDisplayProperty("C", 5)
	    ,new SublineageDisplayProperty("D", 6)
	    ,new SublineageDisplayProperty("E", 2)
	    ,new SublineageDisplayProperty("MS", 4)
	    ,new SublineageDisplayProperty("P", 3)
	    ,new SublineageDisplayProperty("polar", 7)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("", 2)
        ,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
    	,new SublineageDisplayProperty("", 2)
	    ,new SublineageDisplayProperty("other", 2)
	    ,new SublineageDisplayProperty("background", 1)
    };

    public final static int 
        MIN_RED = 25000
        , MAX_RED = 100000
        , SPECIAL = 1
        , MAX_TAIL_TIME_PTS = 20
        , MIN_TAIL_TIME_PTS = 0
        , INIT_TAIL_TIME_PTS = 10
        , DEFAULT_TAIL_OPACITY = 50
        , DEFAULT_OVERLAY_MAX_X = 512
        , DEFAULT_OVERLAY_MIN_X = 0 
        , DEFAULT_OVERLAY_MAX_Y = 512 
        , DEFAULT_OVERLAY_MIN_Y = 0 
        , DEFAULT_OVERLAY_MAX_Z = 30 // This is not being used, the default slider position is set in Image3D2 dynamically
        , DEFAULT_OVERLAY_MIN_Z = 0 
        , DEFAULT_OVERLAY_SUBSAMPLE = 2;
    ;

    public final static Color DEFAULT_CUSTOM_TAIL_COLOR = Color.MAGENTA;
   
    //////////////////// END OF DEFAULTS ////////////////////    

    private static Image3DViewConfig instance = null;
    private     SublineageDisplayProperty[]     dispProps;
    private     boolean                         useExpression;
    private     boolean                         showNonExpressing;
    private     boolean                         useExpressionColors;
    private     int                             minRed;
    private     int                             maxRed;

    private     boolean                         showTails;
    private     Color                           customTailColor;
    private     int                             tailOpacity;
    private     int                             tailTimePts;

    private     boolean                         sistersVisible;
    private     String                          currentRotDir; //
    private     int                             lineageCount;

    private     String                          title;

    private     int                             overlayMaxX;
    private     int                             overlayMinX;
    private     int                             overlayMaxY;
    private     int                             overlayMinY;
    private     int                             overlayMaxZ;
    private     int                             overlayMinZ;
    private     int                             overlaySubsample;
    private     boolean                         overlayXYZChanged;
    private     boolean                         showOverlay;
    private     boolean                         useOverlayAutoROI;
    private     boolean                         useOverlayRedChannel;
    private     boolean                         useOverlayGreenChannel;
    private     boolean                         useOverlayBlueChannel;
    private     boolean                         changeOverlayChannel;

    public Image3DViewConfig() {
        //geometryReset();
        propertiesReset();
    }

    public static Image3DViewConfig getInstance() {
        if(instance == null)
            instance = new Image3DViewConfig();
        return instance;
    }

    // set all attributes to their default values
    public void propertiesReset() {
        dispProps = DISP_PROPS;

        useExpression = false;
        showNonExpressing = false;
        showTails = false;
        minRed = MIN_RED;
        maxRed = MAX_RED;

        tailTimePts = INIT_TAIL_TIME_PTS;
        tailOpacity = DEFAULT_TAIL_OPACITY;
        customTailColor = DEFAULT_CUSTOM_TAIL_COLOR;

        currentRotDir = ".";

        resetOverlayDefaults();
        showOverlay = false;
    }
    
    // Returns default array of sublineage dislpay properties
    public static SublineageDisplayProperty[] getDefaultDispProp() {
    	return DEFAULT_DISP_PROPS;
    }

    public void resetOverlayDefaults() {
        overlayMaxX = DEFAULT_OVERLAY_MAX_X;
        overlayMinX = DEFAULT_OVERLAY_MIN_X;
        overlayMaxY = DEFAULT_OVERLAY_MAX_Y;
        overlayMinY = DEFAULT_OVERLAY_MIN_Y;
        overlayMaxZ = DEFAULT_OVERLAY_MAX_Z;
        overlayMinZ = DEFAULT_OVERLAY_MIN_Z;
        overlayXYZChanged = true;
        overlaySubsample = DEFAULT_OVERLAY_SUBSAMPLE;
        useOverlayAutoROI = false;
        showOverlay = true;
        useOverlayRedChannel = true;
        useOverlayGreenChannel = true;
        useOverlayBlueChannel = true;
        changeOverlayChannel = false;
    }

    // returns number of display properties
    public int getNumDispProps() { 
    	return dispProps.length; 
	}

    //////////////////// ACCESSORS //////////////////// 

    // gets SublineageDisplayProperty at specified index in dispProps
    public SublineageDisplayProperty getDispProp(int i) { 
        return this.dispProps[i];
    }
    
    // Returns the index of the "other" field
    public int getOtherIndex() {
    	for (int i = 0; i < getNumDispProps(); i++) {
    		String name = dispProps[i].getName();
    		name = name.toLowerCase();
    		if (name.equals("other"))
    			return i;
    	}
    	return -1;    	
    }

    public SublineageDisplayProperty[] getDispProps() { 
        return this.dispProps;
    }

    public int getMinRed() { return this.minRed; }
    public int getMaxRed() { return this.maxRed; }
    public boolean isUsingExpression() { return this.useExpression; }
    public boolean isUsingExpressionColors() { 
        return this.useExpressionColors; 
    }
    public boolean isShowingNonExpressing() { return this.showNonExpressing; }
    public boolean areSistersVisible() { return this.sistersVisible; }
    public int getLineageCount() { return this.lineageCount; }
    public String getCurrentRotDir() { return this.currentRotDir; }
    public boolean isShowingTails() { return showTails; }
    public Color getCustomTailColor() { return this.customTailColor; }
    public int getTailOpacity() { return this.tailOpacity; }
    public int getTailTimePoints() { return this.tailTimePts; }
    public Color getSisterColor() { return Color.MAGENTA; }
    public String getTitle() { return this.title; }

    public int getOverlayMaxX() { return this.overlayMaxX; }
    public int getOverlayMinX() { return this.overlayMinX; }
    public int getOverlayMaxY() { return this.overlayMaxY; }
    public int getOverlayMinY() { return this.overlayMinY; }
    public int getOverlayMaxZ() { return this.overlayMaxZ; }
    public int getOverlayMinZ() { return this.overlayMinZ; }
    public int getOverlaySubsample() { return this.overlaySubsample; }

    public boolean isOverlayXYZChanged() { return this.overlayXYZChanged; }
    public boolean isShowingOverlay() { return this.showOverlay; }
    public boolean useOverlayAutoROI() { return this.useOverlayAutoROI; }
    public boolean useOverlayRedChannel() { return this.useOverlayRedChannel; }
    public boolean useOverlayGreenChannel() { return this.useOverlayGreenChannel; }
    public boolean useOverlayBlueChannel() { return this.useOverlayBlueChannel; }
    public boolean changeOverlayChannel() { return this.changeOverlayChannel; }
   
    //////////////////// MUTATORS //////////////////// 

    public void setMinRed(int newMin) { this.minRed = newMin; }
    public void setMaxRed(int newMax) { this.maxRed = newMax; }
    public void setUseExpression(boolean useExpr) { 
        this.useExpression = useExpr; 
    }

    public void setUseExpressionColors(boolean useExprCol) { 
        this.useExpressionColors = useExprCol;
    }

    public void setShowNonExpressing(boolean showNonExpr) { 
        this.showNonExpressing = showNonExpr;
    }

    public void setSistersVisible(boolean sisterVisible) { 
        this.sistersVisible = sisterVisible; 
    }

    private void reportDispProps() {
        for(int i = 0; i < getNumDispProps(); i++) {
            System.out.println("dispProp: " + i + "," + 
                getDispProp(i).getName() + "," + 
                getDispProp(i).getLineageNum()); 
        }
    }

    public void setLineageCount(int count) { this.lineageCount = count; }
    public void setCurrentRotDir(String dir) { this.currentRotDir = dir; } 
    public void setShowingTails(boolean show) { this.showTails = show; }

    public void setCustomTailColor(Color c) {
        if(c != null)
            this.customTailColor = c;
    }

    public void setTailOpacity(int op) { this.tailOpacity = op; }
    public void setTailTimePoints(int tp) { this.tailTimePts = tp; }

    public void setTitle(String t) { this.title = t; }

    public void setOverlayMaxX(int o) { this.overlayMaxX = o; }
    public void setOverlayMaxY(int o) { this.overlayMaxY = o; }
    public void setOverlayMaxZ(int o) { this.overlayMaxZ = o; }
    public void setOverlayMinX(int o) { this.overlayMinX = o; }
    public void setOverlayMinY(int o) { this.overlayMinY = o; }
    public void setOverlayMinZ(int o) { this.overlayMinZ = o; }
    public void setOverlaySubsample(int s) { this.overlaySubsample = s; }

    public void setOverlayXYZChanged(boolean c) { this.overlayXYZChanged = c; }
    public void setShowingOverlay(boolean show) { this.showOverlay = show; }
    public void setUseOverlayAutoROI(boolean use) { this.useOverlayAutoROI = use; }
    public void setUseOverlayRedChannel(boolean use) { this.useOverlayRedChannel = use; }
    public void setUseOverlayGreenChannel(boolean use) { this.useOverlayGreenChannel = use; }
    public void setUseOverlayBlueChannel(boolean use) { this.useOverlayBlueChannel = use; }
    public void setChangeOverlayChannel(boolean use) { this.changeOverlayChannel = use; }
}
