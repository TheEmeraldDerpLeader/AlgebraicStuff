

#include <SymExp.hpp>

#include <QuatRep.hpp>
#include <ArbAlgRep.hpp>

#include <SFML/Graphics.hpp>

#include <vector>
#include <iostream>

const float RAD_TO_DEG = 180.0f/3.1415926535f;

//int CPPBindingTest();

//To do:
//Generate algebra multiplication and addition code
//Look for more 2D algs
//Apply linear exclusion substitution to avoid pseudo-isomorphisms
//Figure out how to get zero divisor algs separately, and find isos between them

//view and compare primary exp group of {(1+i),(1+j)} and {(1+i),(2+j)}

//To Do:
//
//Divide output into grid, find inverse of center of grid, use that as initial approx for reps in that grid block, hopefully should allow for high accuracy, less flickering, and only one computation per object

int CPPBindingTest();

int screenx = 800;
int screeny = 800;
float halfscreenx = screeny/2.0f;
float halfscreeny = screeny/2.0f;

sf::RectangleShape CreateLine(float x1, float y1, float x2, float y2)
{
	x1 = (x1*halfscreeny)+halfscreeny; x2 = (x2*halfscreeny)+halfscreeny; y1 = (-y1*halfscreeny)+halfscreeny; y2 = (-y2*halfscreeny)+halfscreeny;
	float dx = x2-x1; float dy = y2-y1;
	sf::RectangleShape rect;
	float ang = std::atan2(dy, dx);
	float wx = std::cos(ang)*2.0f; float wy = std::sin(ang)*2.0f;
	rect.setRotation(ang*RAD_TO_DEG);
	rect.setPosition(x1+wy, y1-wx);
	rect.setSize(sf::Vector2f(std::sqrt((dx*dx)+(dy*dy)), 4.0f));
	return rect;
}

struct Circle 
{
public:
	QuatRep position = QuatRep();
	QuatRep approx = QuatRep(0,0,0,0);
	float radius = 0.05f;

	Circle() {}
	Circle(QuatRep q) { position = q; approx = position.InvGeoExpRR(); }

	sf::CircleShape ToShape(QuatRep camI)
	{
		//approx = QuatRep(0, 0, 0, 0);
		QuatRep test = camI*position;
		approx = (camI*position).InvGeoExpRR(approx,1);
		QuatRep test2 = QuatGeoExp(approx);
		float x = (approx.i*0.25f*halfscreeny)+halfscreeny; float y = (-approx.j*0.25f*halfscreeny)+halfscreeny;
		float radiusC = radius*halfscreeny;
		sf::CircleShape hold(radiusC); hold.setPosition(x-radiusC, y-radiusC);
		return hold;
	}
};
struct Line
{
public:
	QuatRep p1 = QuatRep();
	QuatRep p2 = QuatRep();
	QuatRep approx1 = QuatRep(0,0,0,0);
	QuatRep approx2 = QuatRep(0,0,0,0);

	Line() {}
	Line(QuatRep q1, QuatRep q2) { p1 = q1; approx1 = p1.InvGeoExpRR(); p2 = q2; approx1 = p2.InvGeoExpRR(); }

	sf::RectangleShape ToShape(QuatRep camI)
	{
		approx1 = (camI*p1).InvGeoExpRR(approx1,1);
		approx2 = (camI*p2).InvGeoExpRR(approx2,1);
		return CreateLine(approx1.i,approx1.j,approx2.i,approx2.j);
	}
};

int main()
{
	std::cout << "test";
	
	SymExp bruh;
	bruh.terms.push_back(Product(1, 1, 2));
	std::cout << bruh.ToString() << '\n';
	
	CPPBindingTest();
	return 0;

	sf::RenderWindow window;
	window.create(sf::VideoMode(screenx, screeny), "Display");
	window.setFramerateLimit(60);


	std::vector<float> testIn; testIn.resize(4);
	std::vector<float> testOut; testOut.resize(4);
	testIn[0] = 1; testIn[2] = 0.1f;
 	QuatGeoExpEval(testIn, testOut);
 	QuatGeoExpGradient(testIn, testOut);
	
	QuatRep test1 = QuatGeoExp(0, 1, 0)*QuatGeoExp(1, 0, 0);
	QuatRep test2 = QuatGeoExp(0, 1, 0);
	for (int i = 0; i < 60; i++)
		test2 *= QuatGeoExp(1/60.0f, 0, 0);

	test1 = test1.InvGeoExpRR();
	test2 = test2.InvGeoExpRR();
	
	test2 = QuatRep();
	for (int i = 0; i < 100; i++)
	{
		for (int i = 0; i < 60; i++)
			test2 *= QuatGeoExp(1/60.0f, 0, 0);
		for (int i = 0; i < 60; i++)
			test2 *= QuatGeoExp(-1/60.0f, 0, 0);
	}

	test2 = test2.InvGeoExpRR();

	//rects.push_back(CreateLine(-0.5, -0.5, 0.5, 0));
	//rects.back().setFillColor(sf::Color::White);

	std::vector<sf::RectangleShape> rects;
	
	std::vector<Circle> circles;
	std::vector<Line> lines;
	circles.push_back(Circle(QuatGeoExp(0, 0, 0)));
	circles.push_back(Circle(QuatGeoExp(0.5f, 0, 0)));
	circles.push_back(Circle(QuatGeoExp(0, 0.5f, 0)));
	//circles.push_back(Circle(QuatGeoExp(0.5f, 0, 0)*QuatGeoExp(0, 1.0f, 0)));

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			//circles.push_back(Circle(QuatGeoExp((j/2.0f)-1.0f, (i/2.0f)-1.0f, 0)));
			//circles.back().radius *= 0.25f;
		}
	}

	QuatRep last;
	for (int t = 1; t < 150; t++)
	{
		QuatRep quat = QuatExp(t*0.02, t*3.14159*2.0f/50.0f, 0, 0);
		//std::cout << quat.ToString() << '\n';
		//rects.push_back(CreateLine(last.s/4.0f,last.i/4.0f,quat.s/4.0f,quat.i/4.0f));
		last = quat;
	}
	rects.push_back(CreateLine(0,-1,0,1));
	rects.push_back(CreateLine(-1,0,1,0));

	QuatRep camera;
	//camera *= QuatGeoExp(0.5f, 0, 0);

	sf::Clock fpsClock;

	while (window.isOpen())
	{
		sf::Time dif = fpsClock.getElapsedTime();
		float deltaT = dif.asSeconds();
		fpsClock.restart();

		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		QuatRep moveVec(0, 0, 0, 0);

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
			moveVec += QuatRep(0,1.0f*deltaT, 0, 0);
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
			moveVec += QuatRep(0,-1.0f*deltaT, 0, 0);
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
			moveVec += QuatRep(0, 0, 1.0f*deltaT, 0);
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
			moveVec += QuatRep(0, 0, -1.0f*deltaT, 0);
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::E))
			camera *= QuatGeoExp(0, 0, 1.0f*deltaT);
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q))
			camera *= QuatGeoExp(0, 0, -1.0f*deltaT);

		camera *= QuatGeoExp(moveVec);
		//camera = camera.LocalMove(moveVec.i, moveVec.j);
		QuatRep camI = camera.Inverse();
		
		window.clear();

		QuatRep testPrev;
		for (int i = 0; i < 40; i++)
		{
			QuatRep next = testPrev*QuatGeoExp(0.2f, 0, 0);
			window.draw(CreateLine(testPrev.s/4.0f, testPrev.i/4.0f, next.s/4.0f, next.i/4.0f));
			testPrev = next;
		}

		for (int i = 0; i < rects.size(); i++)
			window.draw(rects[i]);
		for (int i = 0; i < circles.size(); i++)
			window.draw(circles[i].ToShape(camI));
		window.display();
	}

	return 0;
}