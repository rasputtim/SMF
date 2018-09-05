
#ifndef _TUT_GAME_H_
#define _TUT_GAME_H_

#include "Gamecore/game/SGF_Control.h"
#include "Gamecore/game/SGF_Theater.h"
#include "SGF_Global.h"

class CBaseObject;
class SGF::CBitmap;

class CTutGameControl : public SGF::CGameControl
{


public:
//Contructor Destructor
	//! CMenu_v2 actions

    CTutGameControl();
    ~CTutGameControl();
	// Métodos virtuais herdados da classe Base
	int runMenu() {return true;};
	bool loadConfiguration(){return true;};  //De fato as configurações devems er carregadas antes de inicializar a classe

	//Métodos locais
	void initWindow();
	void setQuit( bool q ) { quit = q; }
	inline bool getQuit() const { return quit; }
	SGF::CBitmap* getWork() { return work;}
	inline SGF::CBitmap *getScreen() {return work;}

protected:


private:
	SGF::CInputManager *pCurInputManager;  //! inicializa Inpout Manager
    //////////////
    //Global
    //////////////
	bool  quit;       //if true then quit application
	SGF::CBitmap *work;      //screen bitmap

};

class Game :public SGF::CTheater
{
public:
    Game();
    ~Game();
	enum Actions{
			Up,
			Down,
			Left,
			Right,
			Select,
			Back,
			Cancel,
			Modify,
			Pause,
			ChangeSkin,
			InitEditor,
			SaveEditor,
			ShowFPS,
			ActionN,
			ChangeLevel,
			ActionH,
			Click,
			Help,
			TogleSound
		};

	virtual void act() ;
	virtual void draw( SGF::CBitmap * work=NULL );
	virtual bool isFinished() const {return getGameControl()->getQuit();};

	virtual void reloadLevel(){}; //evitar reclamação do compilador


	/* upper left hand corner of the screen */
	virtual int getX() {return 0;};
	virtual int getY() {return 0;};
	 /* set to paused */
	virtual void pause() {}; //para evitar do compilador reclamar a falta de um método virtual

	virtual bool isPaused() { return false;}; //para evitar do compilador reclamar a falta de um método virtual


	void setGameControl(CTutGameControl * Game) { gameControl=Game;}
	inline CTutGameControl *getGameControl() const { return gameControl;}
	void addMessage(SGF::Network::Message m, SGF::Network::Socket from = 0, SGF::Network::Socket to = 0){};


    void processInput(int k, int ix = -1, int iy = -1);
    void processLogic();
    void logicGame();

    void gameInit();

    void toggleFps() { show_fps = !show_fps; }

    void clearHscore();

    bool isinit;

    //! Game Loop Process everything
	static void run(Game & thisGame);
	void doAction(Actions Action, int mousex=0, int mousey=0);

private:
	CTutGameControl * gameControl;
    int
            state,
            counter;  //incrementado a cada vez que executa o método render

    std::string fps;

   unsigned int ticks;


    ///////////////////////
    // GAME OBJECTS
    //////////////////////

    int
            key,
            floatingX,
            floatingY,
            namecol[3];

    unsigned int
            time,
            oldtime,
            ghosttick,
            pausetick;

    bool
            inputwaiting,
            gamestarted,
            renderisbusy;
	static bool show_fps;
    std::string
            num[10],
            name;

    //////////////////////////////////
    // EDITOR OBJECTS
    //////////////////////////////////

    int
            activetool,
            mouseX,
            mouseY;

};

#endif
