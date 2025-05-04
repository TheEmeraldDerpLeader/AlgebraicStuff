

#include <SymExp.hpp>

#include <QuatRep.hpp>
#include <ArbAlgRep.hpp>
#include <ApproxGrid.hpp>

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
//Test random 3D QuatReps to see if reasonable results can be found
//If that works, try 4D

int CPPBindingTest();

int screenx = 800;
int screeny = 800;
float halfscreenx = screeny/2.0f;
float halfscreeny = screeny/2.0f;

//ApproxGrid grid;'
LoopApprox loopApprox;

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
	sf::Color color = sf::Color::White;

	Circle() {}
	Circle(QuatRep q) { position = q; approx = position.InvGeoExp(); }

	sf::CircleShape ToShape(QuatRep camI, QuatRep camR)
	{
		//approx = QuatRep(0, 0, 0, 0);
		QuatRep test = camI*position;
		//QuatRep testApprox = grid.GetApprox(test);
		QuatRep testApprox = loopApprox.TryLoopQuot(approx);
		approx = (camI*position).InvGeoExp(testApprox,1);
		QuatRep test2 = QuatGeoExp(approx);
		QuatRep approxAdj = approx.MoveAdjust(camR);
		float x = (approxAdj.i*0.25f*halfscreeny)+halfscreeny; float y = (-approxAdj.j*0.25f*halfscreeny)+halfscreeny;
		float radiusC = radius*halfscreeny;
		sf::CircleShape hold(radiusC); hold.setPosition(x-radiusC, y-radiusC);
		hold.setFillColor(color);
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
	Line(QuatRep q1, QuatRep q2) { p1 = q1; approx1 = p1.InvGeoExp(); p2 = q2; approx1 = p2.InvGeoExp(); }

	sf::RectangleShape ToShape(QuatRep camI)
	{
		approx1 = (camI*p1).InvGeoExp(approx1,1);
		approx2 = (camI*p2).InvGeoExp(approx2,1);
		return CreateLine(approx1.i,approx1.j,approx2.i,approx2.j);
	}
};

int main()
{
	std::cout << "test";
	
	//CPPBindingTest();
	//return 0;
	
	loopApprox.Generate();

	QuatRep testyyyy = loopApprox.TryLoopQuot(loopApprox._movBase[5]*glm::sqrt(1.0f/loopApprox._movBase[5].SqrMag())*loopApprox.loopAssocQuotMag[5]*1.2f);

	SymExp bruh;
	bruh.terms.push_back(Product(1, 1, 2));
	std::cout << bruh.ToString() << '\n';
	

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

	test1 = test1.InvGeoExp();
	test2 = test2.InvGeoExp();
	
	test2 = QuatRep();
	for (int i = 0; i < 100; i++)
	{
		for (int i = 0; i < 60; i++)
			test2 *= QuatGeoExp(1/60.0f, 0, 0);
		for (int i = 0; i < 60; i++)
			test2 *= QuatGeoExp(-1/60.0f, 0, 0);
	}

	test2 = test2.InvGeoExp();

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

	QuatRep last(1,0,0,0);
	for (int t = 1; t < 150; t++)
	{
		QuatRep quat = QuatGeoExp(t*0.02, 0, 0);
		//std::cout << quat.ToString() << '\n';
		rects.push_back(CreateLine(last.s/4.0f,last.i/4.0f,quat.s/4.0f,quat.i/4.0f));
		last = quat;
	}
	for (int i = 0; i < loopApprox.loopAssocQuotMag.size() - 2; i++)
	{
		if (i == 48 || i == 49)
			continue;
		if (loopApprox.loopAssocQuotMag[i] > 0 && loopApprox.loopAssocQuotMag[i+2] > 0)
		{
			glm::vec2 p1 = glm::normalize(glm::vec2(loopApprox._movBase[i].i,loopApprox._movBase[i].j)) * loopApprox.loopAssocQuotMag[i]*1.1f;
			glm::vec2 p2 = glm::normalize(glm::vec2(loopApprox._movBase[i+2].i,loopApprox._movBase[i+2].j)) * loopApprox.loopAssocQuotMag[i+2]*1.1f;
			rects.push_back(CreateLine(p1.x/4.0f,p1.y/4.0f,p2.x/4.0f,p2.y/4.0f));
			rects.back().setFillColor(sf::Color(0,128,0));

			p1 = glm::vec2(loopApprox.loopAssocQuot[i].i, loopApprox.loopAssocQuot[i].j);
			p2 = glm::vec2(loopApprox.loopAssocQuot[i+2].i, loopApprox.loopAssocQuot[i+2].j);
			rects.push_back(CreateLine(p1.x/4.0f,p1.y/4.0f,p2.x/4.0f,p2.y/4.0f));
			rects.back().setFillColor(sf::Color(128,0,128));
		}
	}
	{
		for (int i = 0; i < -1; i++)
		{
			int back = loopApprox.loopAssocQuotMag.size()-2;
			if (loopApprox.loopAssocQuotMag[0+i] > 0 && loopApprox.loopAssocQuotMag[back+i] > 0)
			{
				glm::vec2 p1 = glm::normalize(glm::vec2(loopApprox._movBase[0+i].i, loopApprox._movBase[0+i].j)) * loopApprox.loopAssocQuotMag[0+i]*1.1f;
				glm::vec2 p2 = glm::normalize(glm::vec2(loopApprox._movBase[back+i].i, loopApprox._movBase[back+i].j)) * loopApprox.loopAssocQuotMag[back+i]*1.1f;
				rects.push_back(CreateLine(p1.x/4.0f, p1.y/4.0f, p2.x/4.0f, p2.y/4.0f));
				rects.back().setFillColor(sf::Color(0, 128, 0));

				p1 = glm::vec2(loopApprox.loopAssocQuot[0+i].i, loopApprox.loopAssocQuot[0+i].j);
				p2 = glm::vec2(loopApprox.loopAssocQuot[back+i].i, loopApprox.loopAssocQuot[back+i].j);
				rects.push_back(CreateLine(p1.x/4.0f,p1.y/4.0f,p2.x/4.0f,p2.y/4.0f));
				rects.back().setFillColor(sf::Color(128,0,128));
			}
		}
	}
	rects.push_back(CreateLine(0,-1,0,1));
	rects.push_back(CreateLine(-1,0,1,0));

	//grid.Generate(11, 0.125f,23,0.25f);

	QuatRep camera;
	//camera *= QuatGeoExp(0.5f, 0, 0);

	sf::Clock fpsClock;

	float placementTimer = 5;
	float rotateTimer = 0;

	QuatRep debug1 = QuatGeoExp(0.01, 0, 0)*QuatGeoExp(0, 0.001, 0); debug1 = debug1.InvGeoExp();
	QuatRep debug2 = QuatGeoExp(0.01, 0.001, 0); debug2 = debug2.InvGeoExp();
	debug1 = (debug1-debug2)/(0.01*0.001f*2.0f); std::cout << debug1.ToString() << '\n';

	//geodesic generator
	for (float t = 0; t < 20; t++)
	{
		glm::vec2 p = glm::vec2(glm::cos(-3.14159f*2.0f*t/20.0f), glm::sin(-3.14159f*2.0f*t/20.0f))*0.1f;
		QuatRep pRep = QuatGeoExp(p.x, p.y,0);
		for (int i = 0; i < 10; i++)
		{
			circles.push_back(Circle(pRep));
			circles.back().radius *= 0.15f;
			circles.back().color = sf::Color(255*(i/9.0f), 0, 255*((9-i)/9.0f));

			glm::vec2 bestVec(0, 0);
			for (int a = 0; a < 20; a++)
			{
				QuatRep test = pRep*QuatGeoExp(glm::cos(3.14159f*2.0f*a/20.0f)*0.1f, glm::sin(3.14159f*2.0f*a/20.0f)*0.1f, 0);
				QuatRep invTest = test.InvGeoExp();
				if ((invTest.i*invTest.i)+(invTest.j*invTest.j) >(bestVec.x*bestVec.x)+(bestVec.y*bestVec.y))
					bestVec = glm::vec2(invTest.i, invTest.j);
			}
			p = bestVec;
			pRep = QuatGeoExp(p.x, p.y, 0);
		}

		circles.push_back(Circle(pRep));
		circles.back().radius *= 0.25f;
		circles.back().color = sf::Color(255, 0, 0);
	}

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

		placementTimer -= deltaT;
		if (placementTimer <= 0 && circles.size() < 20 && false)
		{
			placementTimer = 0.2f;
			circles.push_back(camera);
			circles.back().radius = 0.025f;
			circles.back().color = sf::Color(127, 127, 0);
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

		//apparently exp(i)exp(j) isn't necessarily exp(i+j) if they commute. Also for Q and E above, neither of them lie on (i+j) line, which Q should since moving in that direction leads to Q. Inv exp should return (-2,-2,0) which is weird 

		//camera *= QuatGeoExp(moveVec);

		if (std::isnan(camera.s))
			camera = QuatRep(1, 0, 0, 0);

		camera = camera.LocalMove(moveVec.i, moveVec.j);
		//camera *= QuatGeoExp(moveVec.i, moveVec.j, 0);
		QuatRep camR = QuatRotExp(camera.InvGeoExp().k);
		QuatRep camI = camera.Inverse();
		
		window.clear();

		rotateTimer += deltaT*0.5f;
		if (rotateTimer >= 3.14159f*2.0f)
			rotateTimer -= 3.14159f*2.0f;

		for (int i = 0; i < rects.size(); i++)
			window.draw(rects[i]);
		for (int i = 0; i < circles.size(); i++)
			window.draw(circles[i].ToShape(camI,camR));
		window.display();
	}

	return 0;
}

//test surface render
/*
std::vector<QuatRep> testCircles;
for (float x = -4; x <= 4.0f; x += 0.2f)
{
	for (float y = -4; y <= 4.0f; y += 0.2f)
	{
		testCircles.push_back(QuatMovExp(x, y));
	}
}
-----
for (int i = 0; i < testCircles.size(); i++)
{
QuatRep p1 = testCircles[i];
p1 *= 0.25f;
QuatRep p = QuatRep(std::cos(rotateTimer)*p1.s - std::sin(rotateTimer)*p1.i, std::cos(rotateTimer)*p1.i + std::sin(rotateTimer)*p1.s, p1.j, p1.k);

if (p.s < -2.0f)
p.s = -2.0f;
if (p.s > 2.0f)
p.s = 2.0f;

sf::CircleShape hold; hold.setPosition(((p.i-1.0f)*halfscreeny)+screenx,((p.j-1.0f)*halfscreeny)+screeny); //use depth for color
hold.setFillColor(sf::Color(255.0f-(255.0f*(p.s+2.0f)/4.0f),255.0f*(p.s+2.0f)/4.0f,255));

hold.setRadius(2.5f);
window.draw(hold);
}
*/