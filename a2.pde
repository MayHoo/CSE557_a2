Table table;
float[] numbers;
float[] angles;
String[] costs;
int margin = 50;
int y_ticks = 7;
String str = "bar";
String current = "bar";
String str1 = "bar";
float ypercent;
float xpercent;
float percent;
int x = 0;
int y = 0;
int xnew = 0;
int ynew = 0;
int xc = 0;
int yc = 0;
int total = 0;
color[] colors = {color(255,0,0), color(0,255,0),color(0,0,255),color(255,255,0),color(255,0,255),color(0, 255, 255), color(160,255,150), color(97,119,171), color(90,103,39),color(125,0,0), color(0,125,0),color(125,125,125), color(255,125,0), color(255,0,125), color(125,0,255), color(255,125,255), color(125, 255, 255), color(255,255,0), color(255,0,255), color(0, 255, 255)};

float a, b, c, d;
float an, bn, cn, dn;
float step = 250;
int frame = 0;
float[] am;
float[] bm;

void setup(){
  //frameRate(1);
  size(800, 640);
  surface.setResizable(true);
  table = loadTable("data.csv", "header"); // data_a2.csv
  println(table.getRowCount() + " total rows in table"); 
  int n = table.getRowCount();
  numbers = new float[n];
  am = new float[n];
  bm = new float[n];
  costs = new String[n];
  angles = new float[n+1];
  angles[0] = 0;
  angles[n] = 2*PI;
  int i = 0;
  for (TableRow row : table.rows()) {
    String cost = row.getString("cost");
    float time = row.getFloat("time");    
    print(i + " ");
    println(cost + ": " + time);
    costs[i] = cost;
    numbers[i] = time;
    total = (int)(total + time);
    //println (total);
    i++;
  }
  for (int j = 1; j<n; j++) {
     angles[j] = angles[j-1] + radians((numbers[j-1]/total)*360);
     //println(angles[j]);
  }//int total = sum(numbers);
}  

void draw(){
    clear();
    background(0);
    xpercent = width/(float)800;
    ypercent = height/(float)640;
    percent = min(xpercent, ypercent);
    fill(255);
    textSize(percent*20);
    textAlign(CENTER, CENTER);
    //current state
    text(str1, xpercent*100, ypercent*23);
    text("cost/occurrence", xpercent*320, ypercent*23);
    //buttons
    rect(xpercent*560, ypercent*10, xpercent*50, ypercent*30);
    rect(xpercent*630, ypercent*10, xpercent*50, ypercent*30);
    rect(xpercent*700, ypercent*10, xpercent*50, ypercent*30);
    fill(0);
    text("bar",xpercent*585,ypercent*23);  
    text("line",xpercent*655,ypercent*23);  
    text("pie",xpercent*725,ypercent*23);  

    float max_number = max(numbers);
    float h = (640-2*margin)/(max_number);
    float y_tick_h = (int) (h*max_number/(y_ticks-1));
    //println(h);
    float y_tick_val = max_number/(y_ticks-1);
    float w = ((800-2*margin)/numbers.length);
    textAlign(RIGHT, CENTER);
    fill(255);
    if(str == str1) frame = 0;
    if(str == "line" || str == "bar") {
      //ticks and axis for line and bar chart
      for(int i=0;i<y_ticks;i++){
        stroke(255,255,255);
        line(xpercent*margin, ypercent*(margin+i*y_tick_h), xpercent*750, ypercent*(margin+i*y_tick_h));
        text((int)(max_number-(i*y_tick_val)), xpercent*30, ypercent*(50+i*y_tick_h));
        stroke(0);
      }
      textAlign(CENTER);
      textSize(12);
      for(int i=0; i<numbers.length; i++){
        text(costs[i], xpercent*(margin+(i+0.5)*w), ypercent*(640-0.5*margin));
      }
      //int[] numbers = { 2, 5, 3, 1, 6, 5, 9, 4, 7, 3, 2, 5, 1, 4, 2, 5 };
    }
    fill(255,255,255);
    for(int i=0;i<numbers.length;i++){
      if(str1 == "bar") {
        fill(255);
        stroke(0);
        rect(xpercent*(i*w+margin), ypercent*(640-margin), xpercent*w, ypercent*(-h*numbers[i]));
        line(xpercent*margin, ypercent*(640-margin), xpercent*(800-margin), ypercent*(640-margin));
        //stroke(0);
        if ((mouseX>xpercent*(i*w+margin))&&(mouseX<xpercent*(i*w+margin+w))&&(mouseY>ypercent*(640-margin-h*numbers[i]))&&(mouseY<ypercent*(640-margin))) {
        //if ((mouseX>i*w+margin)&(mouseX<i*w+margin+w)) {
          fill(colors[i]);
          rect(xpercent*(i*w+margin), ypercent*(640-margin), xpercent*w, ypercent*(-h*numbers[i]));
          fill(255);
          textSize(12);
          String display = "(" + costs[i] + ", " + (int)numbers[i] + ")";
          text(display, xpercent*(i*w+margin+0.5*w), ypercent*(640-margin-h*numbers[i] - 15));
        } 
      } 
      if(str1 == "line") {
        //fill(255);
        stroke(255);
        fill(255);
        ellipse(xpercent*(i*w+margin+0.5*w), ypercent*(640-h*numbers[i]-margin), percent*(20), percent*(20));
        if(i != numbers.length-1) {
          x = (int)(xpercent*(i*w+margin+0.5*w));
          y = (int)(ypercent*(640-h*numbers[i]-margin));
          xnew = (int)(xpercent*((i+1)*w+margin+0.5*w));
          ynew = (int)(ypercent*(640-h*numbers[i+1]-margin));
          stroke(255);
          line(x, y, xnew, ynew);
        }
        if ((mouseX>xpercent*(i*w+margin))&&(mouseX<xpercent*(i*w+margin+w))&&(mouseY>ypercent*(640-margin-h*numbers[i]-10))&(mouseY<ypercent*(640-margin-h*numbers[i]+10))) {
          //fill(colors[i]);
          //rect(i*w+margin, 640-margin, w, -h*numbers[i]);
          fill(255);
          textSize(12);
          String display = "(" + costs[i] + "," + (int)numbers[i] + ")";
          text(display, xpercent*(i*w+margin+0.5*w), ypercent*(640-margin-h*numbers[i] - 15));
          fill(colors[i]);
          ellipse(xpercent*(i*w+margin+0.5*w), ypercent*(640-h*numbers[i]-margin), percent*(20), percent*(20));
        } 
      }
    }
    if(str1 == "pie") {
      float d = (640-2*margin);
      text("cost", percent*(800-125), percent*2*margin);
      text("time", percent*(800-50), percent*2*margin);
      for (int i = 0; i < numbers.length; i++) {
        //float gray = map(i, 0, numbers.length, 0, 255);
        //fill(255);
        stroke(255);
        fill(colors[i]);
        arc(percent*(800-150)/2, percent*640/2, percent*d, percent*d, angles[i], angles[i+1]);
        float posy = mouseY-percent*(640/2);
        float posx = mouseX-(percent*(800-150)/2);
        if(posx <0 && posy<0 && sqrt(posx*posx+posy*posy)<=percent*d/2) {
          if((atan(posy/posx)+PI)>angles[i] && (atan(posy/posx)+PI)<angles[i+1]) {
            arc(percent*(800-150)/2, percent*640/2, percent*1.1*d, percent*1.1*d, angles[i], angles[i+1]);
            fill(0);
            textSize(percent*12);
            text(costs[i], mouseX, mouseY);
            text((int)numbers[i], mouseX+30, mouseY);
          }
        }
        if(posx <0 && posy>0 && sqrt(posx*posx+posy*posy)<=percent*d/2) {
          if((PI-atan(posy/(-posx)))>angles[i] && (PI-atan(posy/(-posx)))<angles[i+1]) {
            arc(percent*(800-150)/2, percent*640/2, percent*1.1*d, percent*1.1*d, angles[i], angles[i+1]);
            textSize(percent*12);fill(0);text(costs[i], mouseX, mouseY);
            text((int)numbers[i], mouseX+30, mouseY);
          }
        }
        if(posx>0 && posy>0 && sqrt(posx*posx+posy*posy)<=percent*d/2) {
          if((atan(posy/posx))>angles[i] && atan(posy/posx)<angles[i+1]) {
            arc(percent*(800-150)/2, percent*640/2, percent*1.1*d, percent*1.1*d, angles[i], angles[i+1]);
            textSize(percent*12);fill(0);text(costs[i], mouseX, mouseY);
            text((int)numbers[i], mouseX+30, mouseY);
          }
        }
        if(posx >0 && posy<0 && sqrt(posx*posx+posy*posy)<=percent*d/2) {
          if(2*PI-atan((-posy)/posx)>angles[i] && (2*PI-atan((-posy)/posx))<angles[i+1]) {
            arc(percent*(800-150)/2, percent*640/2, percent*1.1*d, percent*1.1*d, angles[i], angles[i+1]);
            textSize(percent*12);fill(0);text(costs[i], mouseX, mouseY);
            text((int)numbers[i], mouseX+30, mouseY);
          }
        }
        fill(colors[i]);textSize(percent*20);
        text(costs[i], percent*(800-125), percent*(3*margin+i*30));
        text((int)numbers[i], percent*(800-75), percent*(3*margin+i*30));
      }
    }
    
    if(str1 != "bar" && str1 != "line" && str1 != "pie") {
      if(current == "bar" && str == "line") {
        for (int j = 0; j<numbers.length; j++) {
          a = j*w+margin;
          an = (j+0.5)*w+margin;
          c = w;
          cn = 10;
          d = h*numbers[j];
          dn = 10;
          if(a<an && c>cn && d>dn) {
            a = a + frame*(an-a)/step;
            c = c - frame*(c-cn)/step;
            d = d - (frame)*(d-dn)/step;
            if (d<10) d = 10;
            if (c<10) c = 10;
          }
          fill(255);
          transform(current, str, j, a, c, d, frame);
        }
      }
      if(current == "line" && str == "bar") {
        for (int j = 0; j<numbers.length; j++) {
          an = j*w+margin;
          a = (j+0.5)*w+margin;
          cn = w;
          c = 10;
          dn = h*numbers[j];
          d = 10;
          if(a>an && c<cn && d<dn) {
            a = a - frame*(a-an)/step;
            c = c + frame*(cn-c)/step;
            d = d + (frame)*(dn-d)/step;
          }
          fill(255);
          transform(current, str, j, a, c, d, frame);
        }
      }
      if ((current == "line"||current == "bar") && str == "pie") {
        for (int j = 0; j<numbers.length; j++) {
          a = (j+0.5)*w+margin;
          b = 640-h*numbers[j]-margin;
          c = 10;
          d = 10;
          an = (325)+((640-2*margin)/2)*cos(angles[j]);
          bn = (640/2)+((640-2*margin)/2)*sin(angles[j]);
          println(an);
          //if(frame<step-20){
          //if(a>(800-150-640+2*margin)/2 && a<(325+270) && b>50 && b<590) {
            a = a+frame*(an-a)/(step-20);
            b = b+frame*(bn-b)/(step-20);
          //}
          transform(current, str, j, a, b, d, frame);
        }
      }
      if (current == "pie" && str != "pie") {
        for (int j = 0; j<numbers.length; j++) {
          a = (j+0.5)*w+margin;
          b = 640-h*numbers[j]-margin;
          c = 10;
          d = 10;
          an = (325)+((640-2*margin)/2)*cos(angles[j]);
          bn = (640/2)+((640-2*margin)/2)*sin(angles[j]);
          println(an);
          //if(frame<step-20){
          //if(a>(800-150-640+2*margin)/2 && a<(325+270) && b>50 && b<590) {
            an = an+frame*(a-an)/(step-20);
            bn = bn+frame*(b-bn)/(step-20);
          //}
          transform(current, str, j, an, bn, d, frame);
        }
      }
      //if(str == "pie") {
        //str1 = "pie";
      //}
      //if(current == "pie" && str1 != "pie") {
        //str1 = str;
      //}
    }
   if(frame<step) frame = frame +1;
}

void transform(String current, String change, int i, float a, float c, float d, int frame){
  int max_number = (int) max(numbers);
  float h = (640-2*margin)/(max_number);
  if (current == "bar"){
    if (change == "line") {
      if (frame < step-10) {
        fill(colors[i]);
        rect(xpercent*a, ypercent*(640-h*numbers[i]-margin), xpercent*c, ypercent*d);
      }
      else if (frame < step) {
        ellipse(xpercent*(a), ypercent*(640-h*numbers[i]-margin), percent*(20), percent*(20));
      }
      else str1 = change;
    }
  }
  if (current == "line"){
    if (change == "bar") {
      if (frame < step) {
        fill(colors[i]);
        //println(a);
        rect(xpercent*a, ypercent*(640-h*numbers[i]-margin), xpercent*c, ypercent*d);
      }
      else str1 = change;
    }
  }
  if (current == "line" || current == "bar"){
    if (change == "pie") {
      if (frame < step-20) {
        fill(colors[i]);
        //println(a);
        rect(xpercent*a, ypercent*b, xpercent*10, ypercent*10);
        am[i] = a;
        bm[i] = b;
      }
      else if(frame<step) {
        fill(colors[i]);
        rect(xpercent*am[i], ypercent*bm[i], xpercent*10, ypercent*10);
        stroke(255);
        line(xpercent*am[i],ypercent*bm[i],xpercent*325,ypercent*320);
      } 
      else str1 = change;
    }
  }
  if (current == "pie"){
    if (change != "pie") {
      if (frame < step-20) {
        fill(colors[i]);
        //println(a);
        rect(xpercent*a, ypercent*c, xpercent*10, ypercent*10);
        am[i] = a;
        bm[i] = b;
      }
      //else if(frame<step) {
      //  fill(colors[i]);
      //  rect(xpercent*am[i], ypercent*bm[i], xpercent*10, ypercent*10);
      //  stroke(255);
      //  line(xpercent*am[i],ypercent*bm[i],xpercent*325,ypercent*320);
      //} 
      else str1 = change;
    }
  }
}

void mouseClicked() {
  current = str;
  if ((mouseX>xpercent*560)&(mouseX<xpercent*610)&(mouseY>ypercent*10)&(mouseY<ypercent*40)) {
    if (str != "bar") {
      str = "bar";        
    } 
  }
  if ((mouseX>xpercent*630)&(mouseX<xpercent*680)&(mouseY>ypercent*10)&(mouseY<ypercent*40)) {
    if (str != "line") {
      str = "line";
    } 
  }
  if ((mouseX>xpercent*700)&(mouseX<xpercent*750)&(mouseY>ypercent*10)&(mouseY<ypercent*40)) {
    if (str != "pie") {
      str = "pie";
    } 
  }
  if (current == str) {
    str1 = str;
  }
  else {
    str1 = current + "->" + str;   
  }
}