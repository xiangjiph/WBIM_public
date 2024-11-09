String serialString = "";
String serialInt = "";

int LCD_ON = 3;
int BRIX_READ = 2;
int RED_LED = 11;
int LCD_STATE = 10;
int VALVE = 9;

int state_Q = 0;

bool LED_Q = false;

//bool on_Q = false;

void setup() {
  Serial.begin(9600);

  pinMode(LCD_ON, OUTPUT);
  pinMode(BRIX_READ, OUTPUT);
  pinMode(RED_LED, OUTPUT);
  pinMode(VALVE, OUTPUT);
 
  pinMode(LCD_STATE, INPUT);

  digitalWrite(LCD_ON, LOW);
  digitalWrite(BRIX_READ, LOW);
  digitalWrite(RED_LED, LOW);
  digitalWrite(VALVE, LOW);

}

void loop() {
  if (Serial.available() > 0) {
    serialString = Serial.readString();
    serialString.trim();

    if(serialString=="LCD"){
      /*if(on_Q){
        Serial.println("Already on!");
      }
      else{*/
        turn_onOff_LCD();
        //Serial.println("Turned on.");
      /*}
      on_Q = true;*/
    }
    /*else if(serialString=="OFF"){
      if(!on_Q){
        Serial.println("Already off!");
      }
      else{
        turn_on_LCD();
        Serial.println("Turned off.");
      }
      on_Q = false;
    }*/
    else if(serialString=="READ"){
      measure_brix();
    }
    else if(serialString=="LED"){
      turn_onOff_LED();
    }
    else if(serialString=="STATE"){
      probe_LCD();
    }
    else if(serialString.substring(0,5)=="VALVE"){
      serialInt = serialString.substring(6);
      open_valve_timed(serialInt.toInt());
    }
    else Serial.println(serialString);
  }
}

void probe_LCD() {
  state_Q = digitalRead(LCD_STATE);
  Serial.println(state_Q);
}

void turn_onOff_LCD() {
  digitalWrite(LCD_ON,HIGH);
  delay(500);
  digitalWrite(LCD_ON,LOW);
  delay(500);
}

void measure_brix(){
  digitalWrite(BRIX_READ, HIGH); 
  delay(500);
  digitalWrite(BRIX_READ, LOW);
  delay(500);
}

void turn_onOff_LED(){
  if(!LED_Q){
    //Serial.println("Turning LED ON");
    analogWrite(RED_LED, 150);
  }
  else{
    //Serial.println("Turning LED OFF");
    digitalWrite(RED_LED, LOW);
  }
  LED_Q = !LED_Q;
}

void open_valve_timed(int ms_time){
  digitalWrite(VALVE, HIGH);
  delay(ms_time);
  digitalWrite(VALVE, LOW);
}