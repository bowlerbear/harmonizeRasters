# Function to read and harmonize global raster layers#
######################################################

#Main arguments is:

#x is a dataframe with 3 columns or list with 3 elements
#myfolder - folder containing original raster layer
#myproj - projection of original raster layer
#name - simple name to call raster layer in output (not necessary....)

#Additional optional arguments:
#newres = desired resolution in units of degrees or metres depending on type of projection
#newextent = desired extent in degrees
#timeperiod = desired time period if raster is a stack
#newproj = desired CRS geographic or projection
#these are examples:
#newproj<-"+proj=eck4 +datum=WGS84"
#newproj<-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
#plot = whether to plot the raster layers as part of the output
#summary = whether to summarise multiple input rasters as the "Mean" or the "Trend"

#default is to harmonize and/all raster layers/bands to a global extent, 1 degree resolution, 
#in geographic projection

#Example use
#whan myrasters dataframe have just one row (i.e. one raster)
outout<-harmonizeRasters(myrasters)
output<-harmonizeRasters(myrasters,newres=500000,newproj="+proj=eck4 +datum=WGS84",plot=TRUE)

#apply to mutliple rasters and combine them as a stack
output<-lapply(1:nrow(myrasters),function(i)harmonizeRasters(myrasters[i,]))
output<-stack(output)

##############################################################################################

harmonizeRasters<-function(x, newres=1,newextent=extent(-180, 180, -90, 90),timeperiod=NULL, 
                           newproj="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs",
                           plot=FALSE,summary=NULL){
  
  #libraries we need
  require(raster)
  require(rgdal)
  require(gdalUtils)
  require(ncdf4)
  
  
  #Extract bits of x for later use in function
  name<-as.character(x["name"])
  myproj<-as.character(x["myproj"])  
  myfolder<-as.character(x["myfolder"]) 
  
  #Create reference raster layer that raster will be reprojected onto
  
  #set up a normal 1 degree degree grid
  ref<-newextent
  ref<-raster(ref)
  res(ref)<-1
  values(ref)<-1#dummy value
  projection(ref)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
  
  #project into whatever is the desired end projection resolution
  refProj<-projectRaster(ref, crs=newproj, res=newres)
  #set dummy values again???
  #values(ref)<-1#dummy value
  refProj<-trim(refProj)  
  
  #Read in original raster layer
  
  #get file types in directory
  files<-list.files(myfolder)
  filetypes<-as.character(sapply(files,function(x)substr(x,nchar(x)-2,nchar(x))))
  
  #if file type is asci or tif
  if(any(filetypes%in%c("tif","asc","wcs",".nc"))){
    myfile<-files[which(filetypes%in%c("tif","asc","wcs",".nc"))]  
    r<-stack(paste(myfolder,myfile,sep="/"))
  }else {
    temp <- new("GDALReadOnlyDataset",myfolder)
    temp2<-asSGDF_GROD(temp)
    r <- stack(temp2)
  }
  
  
  #Look at features of raster - not sure whether we want to print or record these somewhere????
  #r is the original raster
  #str(r)
  #extent(r)
  #crs(r)
  #res(r)
  
  #check whether we need to specific the crs projection of the raster (i..e, if R doesn"t)
  #automatically read it from the files  
  if(is.na(crs(r))){
    projection(r)=myproj
  }
  
  #Dealing with NA values - maybe need to do something here????
  #raster::NAvalue(r) <- -1000 #make sure we don't have any super small NAs
  
  
  #If the original raster has multiple time points, subset to time period of interest
  #I have only tried this with a couple so may not work with others at the moment
  if(!is.null(timeperiod)){
    require(lubridate)
    dates<-as.character(sapply(names(r),function(x)substr(x,2,nchar(x))))
    dates<-gsub("\\.","/",dates)
    dates<-as.Date(dates)
    years=year(dates)
    r <- subset(r, names(r)[years%in%timeperiod]) 
  }
  
  #Again, if the raster has multiple layers, do we just want a summary, i.e., the mean or trend over time
  if(!is.null(summary)){
    if(summary=="Trend"){
      time <- 1:nlayers(r)
      
      fun <- function(x) {
        m <- NA
        try( m <- lm(x ~nlayers)$coefficients[2] ,silent=T)
        m
      }
      
      r <- calc(r, fun)
      
    }else if(summary=="Mean"){
      r <- calc(r, mean)  
    }}
  
  
  #Clip raster to extent of reference grid
  extentGrid <- projectExtent(ref, crs(r))
  r<-crop(r,extentGrid)
  
  
  #Project raster onto the reference grid of raster
  rProj <- projectRaster(r, refProj) 
  
  #Plotting all rasters
  if (plot==TRUE){
    require(ggplot2)
    require(rasterVis)
    
    #original raster
    gplot(r) + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = 'white', high = 'blue') +
      coord_equal()+theme_bw()+ggtitle(name)
    ggsave(paste(name,"png",sep="."))
    
    #reprojected/resampled raster
    gplot(rProj) + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = 'white', high = 'blue') +
      coord_equal()+theme_bw()+ggtitle(paste0(name,"Proj"))
    ggsave(paste(paste0(name,"Proj"),"png",sep="."))
    
    #reprojected/resampled and scaled raster
    #Optional - rescaling raster values between 0 and 1
    #this is just for plotting and the rescaling makes it easier for comparion between diff rasters
    #rProjS<- calc(rProj, function(x)log(x+1))
    rProjS<-rProj
    MIN<-cellStats(rProjS,min)
    MAX<-cellStats(rProjS,max)
    rProjS<- calc(rProjS, function(x)(x-MIN)/(MAX-MIN))
    
    gplot(rProjS) + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = 'white', high = 'blue') +
      coord_equal()+theme_bw()+ggtitle(paste0(name,"ProjS"))
    ggsave(paste(paste0(name,"ProjS"),"png",sep="."))
  }
  
  #Return new layer
  return(rProj)
  
}
