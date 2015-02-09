library("qtbase")

type = NULL
sample = NULL
target = NULL

qsetClass("Window", Qt$QWidget, function(parent = NULL) {
  super(parent)
  
  grid <- Qt$QGridLayout()
  ## NOTE: the layout does not take ownership of the widgets, so we
  ## need to assign the layout to our widget up-front. Our widget then
  ## takes ownership of the widgets in the layout.
  setLayout(grid) 
  #grid$addWidget(createFirstExclusiveGroup(), 0, 0)
  grid$addWidget(createTypeGroup(), 0, 0)
  grid$addWidget(createTargetGroup(), 1, 0)
  grid$addWidget(createSampleGroup(), 2, 0)
  grid$addWidget(createDownloadGroup(), 3, 0)
  #grid$addWidget(createComboBoxGroup(), 0, 1)
  #grid$addWidget(createPushButtonGroup(), 1, 1)
  
  setWindowTitle("Encode downloader")
  resize(200, 320)
})

qsetMethod("createPushButtonGroup", Window, function() {
  groupBox <- Qt$QGroupBox("&Push Buttons")
  groupBox$setCheckable(TRUE)
  groupBox$setChecked(TRUE)
  
  pushButton <- Qt$QPushButton("&Normal Button")
  toggleButton <- Qt$QPushButton("&Toggle Button")
  toggleButton$setCheckable(TRUE)
  toggleButton$setChecked(TRUE)
  flatButton <- Qt$QPushButton("&Flat Button")
  flatButton$setFlat(TRUE)
  
  popupButton <- Qt$QPushButton("Pop&up Button")
  menu <- Qt$QMenu(this)
  menu$addAction("&First Item")
  menu$addAction("&Second Item")
  menu$addAction("&Third Item")
  menu$addAction("F&ourth Item")
  popupButton$setMenu(menu)
  
  newAction <- menu$addAction("Submenu")
  subMenu <- Qt$QMenu("Popup Submenu")
  subMenu$addAction("Item 1")
  subMenu$addAction("Item 2")
  subMenu$addAction("Item 3")
  newAction$setMenu(subMenu)
  
  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(pushButton)
  vbox$addWidget(toggleButton)
  vbox$addWidget(flatButton)
  vbox$addWidget(popupButton)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)
  
  groupBox
}, "private")


qsetMethod("createTypeGroup", Window, function() {
     groupBox <- Qt$QGroupBox("Type")
     groupBox$setCheckable(TRUE)
     groupBox$setChecked(FALSE)
   
     radio1 <- Qt$QRadioButton("Experiment")
     radio2 <- Qt$QRadioButton("Web page")
     radio3 <- Qt$QRadioButton("Publication")
     radio1$setChecked(TRUE)

     vbox <- Qt$QVBoxLayout()
     vbox$addWidget(radio1)
     vbox$addWidget(radio2)
     vbox$addWidget(radio3)
     vbox$addStretch(1)
     groupBox$setLayout(vbox)
   
     groupBox
   }, "private")

qsetMethod("createTargetGroup", Window, function() {
  groupBox <- Qt$QGroupBox("Target")
  groupBox$setCheckable(TRUE)
  groupBox$setChecked(FALSE)
  
  radio1 <- Qt$QRadioButton("histone modification")
  radio2 <- Qt$QRadioButton("transcription factor")
  radio3 <- Qt$QRadioButton("control")
  radio1$setChecked(TRUE)
  
  
  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(radio1)
  vbox$addWidget(radio2)
  vbox$addWidget(radio3)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)
  
  groupBox
}, "private")

setSampleTissue <- function (bool){
  if(bool){
    sample = "tissue"
  }
  print(sample)
}

setSampleCell <- function (bool){
  if(bool){
    sample = "primary cell"
  }
  print(sample)
}

setSample <- function (bool){
  if(!bool){
    sample = NULL
  }
  print(sample)
}

qsetMethod("createSampleGroup", Window, function() {
  groupBox <- Qt$QGroupBox("Sample")
  groupBox$setCheckable(TRUE)
  groupBox$setChecked(FALSE)
  
  radio1 <- Qt$QRadioButton("tissue")
  radio2 <- Qt$QRadioButton("primary cell")
  radio1$setChecked(TRUE)
  qconnect(groupBox, "clicked(bool)", function(bool) setSample(bool) )
  qconnect(groupBox, "clicked(bool)", function(bool) setSampleCell(!bool) )
  qconnect(radio1, "toggled(bool)", function(bool) setSampleTissue(bool) )
  qconnect(radio2, "toggled(bool)", function(bool) setSampleCell(bool) )
  
  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(radio1)
  vbox$addWidget(radio2)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)
  
  groupBox
}, "private")

qsetMethod("createComboBoxGroup", Window, function() {
  groupBox <- Qt$QGroupBox("&ComboBox")
  groupBox$setCheckable(TRUE)
  groupBox$setChecked(TRUE)
  
  cbTypeLabel <- Qt$QLabel("Type")
  cbType <-Qt$QComboBox()
  cbType$addItem("Experiment",1)
  cbType$addItem("Publications",2)
  cbType$addItem("Web page",3)
  #qconnect(cbType,  SIGNAL("activated(int)", function(row) print(row)))
  #cbType$show() 
  
  cbTargetLabel <- Qt$QLabel("Target")
  cbTarget <-Qt$QComboBox()
  cbTarget$addItem("histone modification",1)
  cbTarget$addItem("transcription factor",2)
  cbTarget$addItem("control",3)
  
  
  cbSampleLabel <- Qt$QLabel("Biosample type")
  cbSample <-Qt$QComboBox()
  cbSample$addItem("tissue",1)
  cbSample$addItem("primary cell",2)
  
  grid <- Qt$QGridLayout()
  grid$addWidget(cbTypeLabel, 0,0)
  grid$addWidget(cbType, 0,1)
  grid$addWidget(cbTargetLabel, 1,0)
  grid$addWidget(cbTarget, 1,1)
  grid$addWidget(cbSampleLabel, 2,0)
  grid$addWidget(cbSample, 2,1)
  
  #grid$addStretch(1)
  
  groupBox$setLayout(grid) 
  groupBox
}, "private")

qsetMethod("createDownloadGroup", Window, function() {
  groupBox <- Qt$QGroupBox()

  downloadBt <- Qt$QPushButton("Download!")
  qconnect(downloadBt, "pressed", function() print("Downloading"))
  
  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(downloadBt)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)
  
  groupBox
}, "private")

Window()$show()