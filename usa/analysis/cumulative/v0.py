import numpy as np

from usa.analysis.base import vrbs

class VisualiseTwoTangent(dict):

    def __init__(self):
        """
            Object to store visualisation data related to the two-tangent method.
        """
        self["complete"] = False
        self["mode"] = None
        self["std_value"] = None

  
def two_tangent_standardisation(
        pd_data : np.ndarray,
        pd_low_threshold : np.float,
        visualise = False, 
        verbose = False
        ):
    """
        Refactor of original standardisation method for modularity and visualisation; the algorithm itself is essentially unchanged.

        Standardisation is a two-staged process.

        Parameters:
        pd_data         : np.ndaray (flat)
    """

    pd_data = pd_data.astype(int)

    vrbs("* I am in two-tangent method...", verbose)

    vis1, vis2 = VisualiseTwoTangent(), VisualiseTwoTangent()
    
    refImageScalars = np.array(pd_data)

    # Stage 1:
    try:
        refImageScalars = refImageScalars[refImageScalars > pd_low_threshold]
        scalarWMF = int(pd_low_threshold)
        scalarmin = np.int64(0)
        scalarmax = np.int64(np.ceil((np.max(refImageScalars))))
        scalar_range = (scalarmax - scalarmin) +1        
        xaxis = np.zeros(scalar_range) 
        yaxis = np.zeros(scalar_range)
        std_value1, std_value2 = None, None
        step2_calculated = False

        for a in range(1,scalar_range+1):            
            xaxis[a-1] = a      
        #print("Calculating FMBV....")
                    
        for a in range(1,np.size(refImageScalars)):
            yaxis[refImageScalars[a]] = yaxis[refImageScalars[a]] + 1
        
        #produce cumulative summing of these values
        #plot the cumulative pdf of the ROI.
        cumsum = np.cumsum(yaxis)
        stg1_cumsum = cumsum


        vis1["x"] = xaxis
        vis1["y_hist"] = yaxis
        vis1["y_cumsum"] = cumsum

        # print(xaxis,cumsum) # !!
        (ar,br) = np.polyfit(xaxis, cumsum, 1)
        linfit = np.polyval([ar,br],xaxis)
        stg1_linfit = linfit

        vis1["y_linfit_0"] = linfit[0], # line of best fit to calculate intersection points
        # fig_data_std_1["y_linfit_0"] = fig_data_std_1["y_linfit_0"][0] # ...[0] to look inside the tuple structure. I'm not sure why we have to do it like this but its the only way I could get it to work???

        #start from where gradient begins
        st = scalarmin               
        #absscia turning pt calculation
        absscia = cumsum - linfit
        stg1_absscia = absscia

        st_i = 0
        while cumsum[st_i] < 1:
            st_i+=1
        
        #find tangent points
        absscianum = np.argmax(absscia[scalarWMF:])+scalarWMF 
        pt2 = np.nanmax([st_i,absscianum])
        
        if cumsum[pt2] < linfit[pt2]:
            while cumsum[pt2] < linfit[pt2]:
                pt2 = pt2+1
        else:
            while cumsum[pt2] > linfit[pt2]:
                pt2 = pt2+1

        #print st
        pt1 = np.nanmax([st_i,scalarWMF])
        while cumsum[pt1] < linfit[pt1]:               
            pt1 = pt1+1

        stg1_pt2 = pt2
        stg1_pt1 = pt1

        vis1["pt_a"] =  xaxis[pt1]
        vis1["pt_b"] = xaxis[pt2]

        #1st tangent linear regression on pts +/- 10 away..
        #bounds checking
        if pt1 > scalarmin+11:                
            coeffspt1 = np.polyfit(xaxis[pt1-10:pt1+10], cumsum[pt1-10:pt1+10],1)
        else:
            diff = 10-pt1
            coeffspt1 = np.polyfit(xaxis[pt1-10+diff:pt1+10], cumsum[pt1-10+diff:pt1+10],1)
            
        #run the tangent from the point of intersection up to the 2nd intersection point
        pt1linfit = np.polyval(coeffspt1,xaxis[pt1:pt2] )
        stg1_pt1linfit = pt1linfit

        vis1["y_linfit_a"] = np.polyval(coeffspt1,xaxis)

        if pt2 < scalarmax-11 :                
            coeffspt2 = np.polyfit(xaxis[pt2-10:pt2+10], cumsum[pt2-10:pt2+10],1)
            pt2linfit = np.polyval(coeffspt2,xaxis[pt1:pt2+10] )

        else:
            diff = 10-pt2
            coeffspt2 = np.polyfit(xaxis[pt2-10:pt2+10-diff], cumsum[pt2-10:pt2+10-diff],1)    
            pt2linfit = np.polyval(coeffspt2,xaxis[pt1:pt2+10-diff] )
        stg1_pt2linfit = pt2linfit

        vis1["y_linfit_b"] = np.polyval(coeffspt2,xaxis)

        i=0
        while pt1linfit[i] < pt2linfit[i]:
            i = i+1
        i = pt1 + i

        pos_ndg = 10+scalarWMF
        neg_ndg = 10-scalarWMF
        
        coeffsabs = np.polyfit( xaxis[np.argmax(absscia[scalarWMF:])-neg_ndg:np.argmax(absscia[scalarWMF:])+pos_ndg], absscia[np.argmax(absscia[scalarWMF:])-neg_ndg:np.argmax(absscia[scalarWMF:])+pos_ndg] ,2)
        polyfitabs = np.polyval(coeffsabs,absscia[np.argmax(absscia[scalarWMF:]-neg_ndg):np.argmax(absscia[scalarWMF:])+pos_ndg])  

        std_value1 = min([absscianum,i])  #+ scalarWMF

        vrbs("std_value_1 = "+ str(std_value1), verbose)

        vis1["std_value"] = std_value1
        vis1["complete"] = True
        
    except Exception as e:
        print(e)
        std_value1 = None

    # Stage 2:
    try:
        #convert vector of n-scalars into binned sv->256 length vector
        sv = std_value1        
        yaxis = np.zeros(scalar_range-sv)
        for a in range(1,np.size(refImageScalars)):
            if refImageScalars[a] > sv:
                yaxis[refImageScalars[a] - sv] = yaxis[refImageScalars[a]- sv] +1

        xaxissv = np.zeros(scalar_range-sv) 
        for a in range(1,(scalar_range+1)-sv):            
            xaxissv[a-1] = a 
        stg2_xaxis = xaxissv
        #produce cumulative summing of these values
        cumsum = np.cumsum(yaxis)
        stg2_cumsum = cumsum

        (ar,br) = np.polyfit(xaxissv, cumsum, 1)
        linfit = np.polyval([ar,br],xaxissv)
        stg2_linfit = linfit

        #absscia turning pt calculation
        absscia = cumsum - linfit
        stg2_absscia = absscia
        abssicianum = np.argmax(absscia)

        #find tangent points
        pt2 = abssicianum
        while cumsum[pt2] > linfit[pt2]:                        
            if pt2 >= min(len(cumsum)-5, len(linfit)-5):
                break
            pt2 = pt2+1
            #print(pt2, abssicianum, len(linfit), len(cumsum))
        stg2_pt2 = pt2

        pt1 = 0
        while cumsum[pt1] < linfit[pt1]:            
            pt1 = pt1+1
        stg2_pt1 = pt1

        #out of bounds checks -- lower bounds
        if pt1 > scalarmin+10:                
            coeffspt1 = np.polyfit(xaxis[pt1-10:pt1+10], cumsum[pt1-10:pt1+10],1)
        else:
            diff = 10-pt1
            coeffspt1 = np.polyfit(xaxis[pt1-10+diff:pt1+10], cumsum[pt1-10+diff:pt1+10],1)

        #1st tangent linear regression on pts +/- 10 away..
        #coeffspt1 = polyfit(xaxissv[pt1-10:pt1+10], cumsum[pt1-10:pt1+10],1)
        #run the tangent from the point of intersection up to the 2nd intersection point
        pt1linfit = np.polyval(coeffspt1,xaxissv[pt1:pt2] )
        stg2_pt1linfit = pt1linfit

        if pt2 < scalarmax-11 :                
            coeffspt2 = np.polyfit(xaxissv[pt2-10:pt2+10], cumsum[pt2-10:pt2+10],1)
            pt2linfit = np.polyval(coeffspt2,xaxissv[pt1:pt2+10] )
        else:
            diff = 10-pt2
            coeffspt2 = np.polyfit(xaxissv[pt2-10:pt2+10-diff], cumsum[pt2-10:pt2+10-diff],1)
            pt2linfit = np.polyval(coeffspt2,xaxissv[pt1:pt2+10-diff] )
        stg2_pt2linfit = pt2linfit

        #need to fit 2nd order polynomial to this point +/-
        i=0 
        while pt1linfit[i] < pt2linfit[i]:
            i = i+1
        i = pt1 + i            

        #coeffsabs = np.polyfit( xaxissv[np.argmax(absscia)-10:np.argmax(absscia)+10], absscia[np.argmax(absscia)-10:np.argmax(absscia)+10] ,2)
        #polyfitabs = np.polyval(coeffsabs,absscia[np.argmax(absscia)-10:np.argmax(absscia)+10])
        std_value2 =  np.nanmin([abssicianum,i]) + std_value1
        step2_calculated = True

        vrbs("std_value_2 = " + str(std_value2), verbose)

        vis2 = {
            "complete": True,

            "x": xaxissv,
            "y_hist": yaxis,
            "y_cumsum": cumsum,

            "y_linfit_0": linfit, # line of best fit to calculate intersection points
            "pt_a": xaxissv[pt1],
            "pt_b": xaxissv[pt2],

            "y_linfit_a": np.polyval(coeffspt1,xaxissv),
            "y_linfit_b": np.polyval(coeffspt2,xaxissv),

            "std_value": np.nanmin([abssicianum,i])
        }
        
    except (ValueError, IndexError, TypeError, UnboundLocalError) as e:
        std_value2 = None

    return std_value1, std_value2, vis1, vis2

def normalise(vals, sv):
    vals = vals/float(sv)

    for i in range(np.size(vals)): # clip vals at 1.0
        if vals[i] >= 1.0:
            vals[i] = 1.0

    return np.nanmean(vals)*100 # this is the final fmbv index (as a %)