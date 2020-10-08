# This is just a picture and doesn't depend on data
library(diagram)
#generate screening interval/slip rate diagram
# timeline -> box-box-box-box-> box per stage
# screening events decorate boxes


box_ids<-c("I", "II", "III", "IV")


individual_level<-0.5

time_block<-c(4,2,2,1)
time_block<-time_block/sum(time_block)
total_width<-0.7
local_center<-0.1+total_width*cumsum(c(0,time_block[1:3]))+0.5*time_block
local_radius<-total_width*time_block*0.5
start_box<-local_center-local_radius
end_box<-local_center+local_radius

screen_points<-c(start_box[1]*0.5,local_center[2],(end_box[4]+1)*0.5)

setEPS()
postscript(sprintf("figs/%s_figure_2_screening_interval.eps",date_code) ,
           horizontal = FALSE, onefile = FALSE, paper = "special",height=10,width=10)

openplotmat(main = "Screening interval affects detection rates")

straightarrow(from=c(screen_points[1],0.5),to=c(screen_points[1],0.3),lty=2)
straightarrow(from=c(screen_points[2],0.5),to=c(screen_points[2],0.3),lty=2)
straightarrow(from=c(screen_points[3],0.5),to=c(screen_points[3],0.3),lty=2)

textdiamond(mid=c(screen_points[1],0.3-0.045),radx=0.05,rady=0.045,lab="Screen")
textdiamond(mid=c(screen_points[2],0.3-0.045),radx=0.05,rady=0.045,lab="Screen")
textdiamond(mid=c(screen_points[3],0.3-0.045),radx=0.05,rady=0.045,lab="Screen")

text(0.5,0.75,"Example temporal evolution of one individual",cex=1.5)


for(ii in 1:4){
  textrect(mid=c(local_center[ii],0.5),radx=local_radius[ii],rady=0.045,shadow.size=0,lab=box_ids[ii],cex=1.5)
}

straightarrow(from=c(0,0.5),to=c(start_box[1],0.5))
straightarrow(from=c(end_box[4],0.5),to=c(1.0,0.5))


box()
dev.off()
