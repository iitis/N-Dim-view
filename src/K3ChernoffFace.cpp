#include "K3ChernoffFace.h"

#include "K3Helpers.h"

// visual[0] is Skin_C
void draw_Face_oval(QPainter& painter, int cx, int cy, int cxFace, int cyFace, std::optional<double> feature)
{
	CRGBA mujKolor;
	if (feature.has_value())
		mujKolor = K3_color(feature.value(), 1.0); // ;(0.3, 0.4)
	else
		mujKolor = K3_color(0.0, 0.0);


	int k3R, k3G, k3B, k3A;
	k3R = mujKolor.red(); k3G = mujKolor.green();
	k3B = mujKolor.blue(); k3A = mujKolor.alpha();
	QColor QmujKolor(mujKolor.red(), mujKolor.green(), mujKolor.blue(), mujKolor.alpha());


	QPen pen;
	pen.setWidth(6); //  (k3R + k3G - k3B - k3A);
	pen.setColor(QmujKolor);
	// pen.setColor(QColor(255, 255, 0, 0));
	painter.setPen(pen);
	painter.setBrush(QmujKolor);
	// painter.drawPoint(25, 25);
	// Face oval:
	painter.drawEllipse((cx - cxFace) / 2, (cy - cyFace) / 2, cxFace, cyFace);
}

// visual[1] is Hair_C, visual[7] is Hair_L
QColor draw_Hair(QPainter& painter, int cx, int cy, int cxFace, int cyFace, std::optional<double> feature1, std::optional<double> feature7)
{
	QColor Qhaircolor(qRgba(0,0,0,0));
	if (feature1.has_value() && feature7.has_value()) {
		CRGBA mujKolor = K3_color(feature1.value(), 1.0);
		int k3R, k3G, k3B, k3A;
		k3R = mujKolor.red(); k3G = mujKolor.green();
		k3B = mujKolor.blue(); k3A = mujKolor.alpha();

		QColor Qhaircolor;
		Qhaircolor.setAlpha((int)(mujKolor.alpha()));
		Qhaircolor.setRed((int)(mujKolor.red()));
		Qhaircolor.setGreen((int)(mujKolor.green()));
		Qhaircolor.setBlue((int)(mujKolor.blue()));

		QPen pen;
		pen.setColor(Qhaircolor);
		painter.fillRect((cx - cxFace / 4) / 2, (cy - cyFace) / 2, cxFace / 4, cyFace * (feature7.value()) / 4, Qhaircolor);
	}

	return Qhaircolor;
}

// visual[2] is Eye_S
std::tuple<int, int, int, int, int, int, int, int> draw_Eyes(QPainter& painter, int cx, int cy, int cxFace, int cyFace, std::optional<double> feature)
{
	double factEyes = 0.0;
	if (feature.has_value()) {
		factEyes = (feature.value() + 1.0) * 0.5;
	}

	int yEyes = cy / 2 - cyFace / 3;
	int dxEyes = cxFace / 3;
	int cxEyes = (factEyes * cxFace) / 4;
	int cyEyes = (factEyes * cyFace) / 4;
	int RightExoCanthion, RightEndoCanthion, LeftEndoCanthion, LeftExoCanthion;
	RightExoCanthion = cx / 2 - dxEyes - cxEyes / 2;
	RightEndoCanthion = RightExoCanthion + cxEyes;
	LeftEndoCanthion = cx / 2 + dxEyes - cxEyes / 2;
	LeftExoCanthion = LeftEndoCanthion + cxEyes;

	// pen.setColor(Qt::black);
	painter.setBrush(Qt::white);
	painter.setPen(Qt::black);
	// painter.drawEllipse(cx / 2 - dxEyes - cxEyes / 2, yEyes - cyEyes / 2, cxEyes, cyEyes);
	// painter.drawEllipse(cx / 2 + dxEyes - cxEyes / 2, yEyes - cyEyes / 2, cxEyes, cyEyes);
	painter.drawEllipse(cx / 2 - dxEyes - cxEyes / 2, yEyes - cyEyes / 2, cxEyes, cyEyes);
	painter.drawEllipse(cx / 2 + dxEyes - cxEyes / 2, yEyes - cyEyes / 2, cxEyes, cyEyes);

	return { dxEyes, yEyes, cxEyes, cyEyes, RightExoCanthion, RightEndoCanthion, LeftEndoCanthion, LeftExoCanthion };
}

// visual[3] is Nose_L
void draw_visual_3(QPainter& painter, std::optional<double> feature)
{

}

// visual[4] is Mouth_W
void draw_visual_4(QPainter& painter, std::optional<double> feature)
{

}

// visual[5] is Smile
void draw_visual_5(QPainter& painter, std::optional<double> feature)
{

}

// visual[6] is Frown
void draw_visual_6(QPainter& painter, std::optional<double> feature)
{

}

// visual[7] is Hair_L
void draw_visual_7(QPainter& painter, std::optional<double> feature)
{
	// see: draw_Hair
}

// visual[8] is Face_Elong
std::pair<int, int> get_Face_Elong(QPainter& painter, int cx, int cy, std::optional<double> feature)
{
	int cxFace, cyFace;
	if (feature.has_value()) {
		cxFace = cx * (0.6 - 0.38 * feature.value());
		cyFace = cy * (0.6 + 0.38 * feature.value());
	}
	else {
		cxFace = cx * 0.6;
		cyFace = cy * 0.6;
	}

	return { cxFace, cyFace };
}

// visual[9] is Iris_C
QColor draw_Irises(QPainter& painter, int cx, int cy, int dxEyes, int yEyes, int cxEyes, int cyEyes, std::optional<double> feature)
{
	CRGBA mujKolor;
	if (feature.has_value())
		mujKolor = K3_color(0.5 + 0 * feature.value(), 1.0);
	else
		mujKolor = K3_color(0.0, 0.0);

	int k3R, k3G, k3B, k3A;
	k3R = mujKolor.red(); k3G = mujKolor.green();
	k3B = mujKolor.blue(); k3A = mujKolor.alpha();

	QColor QmujKolor;
	QmujKolor.setAlpha((int)(mujKolor.alpha()));
	QmujKolor.setRed((int)(mujKolor.red()));
	QmujKolor.setGreen((int)(mujKolor.green()));
	QmujKolor.setBlue((int)(mujKolor.blue()));

	painter.setBrush(QmujKolor);
	painter.setPen(QmujKolor);
	painter.drawEllipse(cx / 2 - dxEyes - cxEyes / 4, yEyes - cyEyes / 4, cxEyes / 2, cyEyes / 2);
	painter.drawEllipse(cx / 2 + dxEyes - cxEyes / 4, yEyes - cyEyes / 4, cxEyes / 2, cyEyes / 2);

	return QmujKolor;
}


K3ChernoffFace::K3ChernoffFace(int cx, int cy, std::vector<std::optional<double>> k3params) : CImage(cx, cy, CImage::Format::Format_ARGB32) {
	//  Feature mapping:
		// +0 skin colour (rainbow scale)
		// +1 hair colour (incl. facial hair, if any)
		// +2 - eye size
		// +3 - nose height
		// +4 - mouth width
		// +5 - smile
		// 6 - eybrow frown
		// +7 - hair length
		// +8 - face elongation
		// +9 - iris color

	
	QImage* Qthis = (QImage*)this;
	QPainter painter(Qthis);

	QColor CanvasColor;
	CanvasColor.setAlpha((int)0);
	CanvasColor.setRed((int)255);
	CanvasColor.setGreen((int)255);
	CanvasColor.setBlue((int)111);

	QPen pen;
	pen.setColor(CanvasColor);

	painter.fillRect(0, 0, cx, cy / 2, CanvasColor);

	// get face oval
	auto [cxFace, cyFace] = get_Face_Elong(painter, cx, cy, k3params[8]);
	
	// get face color and draw face oval
	draw_Face_oval(painter, cx, cy, cxFace, cyFace, k3params[0]);

	int yf0 = (cy - cyFace) / 2;
	int yf1 = yf0 + (1.0 * cyFace) / 12;
	int yf2 = yf0 + (2.0 * cyFace) / 12;
	int yf3 = yf0 + (3.0 * cyFace) / 12;
	int yf4 = yf0 + (4.0 * cyFace) / 12;
	int yf8 = yf0 + (8.0 * cyFace) / 12;
	int yf11 = yf0 + (11.0 * cyFace) / 12;


	// Now draw hair:
	QColor Qhaircolor = draw_Hair(painter, cx, cy, cxFace, cyFace, k3params[1], k3params[7]);


	// Now eyes:
	auto [dxEyes, yEyes, cxEyes, cyEyes, RightExoCanthion, RightEndoCanthion, LeftEndoCanthion, LeftExoCanthion] =
		draw_Eyes(painter, cx, cy, cxFace, cyFace, k3params[2]);


	// Irises:
	QColor QmujKolor = draw_Irises(painter, cx, cy, dxEyes, yEyes, cxEyes, cyEyes, k3params[9]);


	// Pupils:
	painter.setBrush(Qt::black);
	painter.setPen(QmujKolor);
	// painter.drawEllipse(cx / 2 - dxEyes - cxEyes / 6, yEyes - cyEyes / 6, cxEyes / 3, cyEyes / 3);
	// painter.drawEllipse(cx / 2 + dxEyes - cxEyes / 6, yEyes - cyEyes / 6, cxEyes / 3, cyEyes / 3);

	if (k3params[0].has_value())
	{ // Digression - eyebrows:
		int yBrowsIn = cy / 2 - cyFace * 2 / 5 - k3params[0].value();
		int yBrowsOut = cy / 2 - cyFace * 2 / 5;
		// pen.setBrush(Qt::blue);
		// pen.setColor(Qt::red);
		pen.setColor(Qhaircolor);
		painter.setPen(pen);
		painter.drawLine(RightExoCanthion, yBrowsOut, RightEndoCanthion, yBrowsIn); //  (int x1, int y1, int x2, int y2);
		painter.drawLine(LeftExoCanthion, yBrowsOut, LeftEndoCanthion, yBrowsIn); //  (int x1, int y1, int x2, int y2);
		// painter.drawLine(10, 10, 300, 300); //  (int x1, int y1, int x2, int y2);
	}


	// Nose:
	if (k3params[3].has_value()) {
		painter.setBrush(Qt::white);
		painter.setPen(Qt::black);
		int cxNose = cxEyes;
		int cyNose = cyFace * (1.0 + k3params[3].value()) / 4;
		painter.drawEllipse((cx - cxNose) / 2, yEyes, cxNose, cyNose);
	}


	// Mouth:
	if (k3params[5].has_value())
	{

		//	int cxMouth = cxFace * (0.4 * (0.5 + k3params[4]));
		//	int cyMouth = cyFace * (0.3 * (abs(k3params[5] - 0.5) + .05));
		//	int cxMouth = cxFace * 0.8;
		//	int cyMouth = cyFace * 0.6;
		pen.setWidth(30);
		pen.setColor(Qt::red);
		painter.setPen(pen);	// painter.drawLine(10, 10, 300, 300);

		//	jeśli rozmiar ust ma być stały i powiazany z rozmiarem buźki, to raczej tak :
			//int cxMouth = cxFace * (0.4 * (0.5 + k3params[4]));
			//int cyMouth = cyFace * (0.3 * (abs(k3params[5] - 0.5) + .05));
		int cxMouth = cxFace / 1.4;
		int cyMouth = cyFace / 3.6;

		// painter.drawLine(10, 10, 300, 300);

		if (k3params[5] > 0.5) {
			painter.drawArc((cx - cxMouth) / 2, (cy - cyMouth) / 2 + cyFace * 0.007, cxMouth, cyMouth, 1, (180 * 16) - 1);
			//painter.drawArc(1, 1, 100, 100, 1, (180 * 16) - 2);
		}
		else
		{
			painter.drawArc((cx - cxMouth) / 2, (cy - cyMouth) / 2 + cyFace * 0.007, cxMouth, cyMouth, (180 * 16) + 2, (180 * 16) - 2);
			//painter.drawArc(1, 1, 100, 100, (180 * 16) + 2, (180 * 16) - 2);
		};


#if 0
		if (k3params[5] > 0.5) {
			painter.drawArc((cx - cxMouth) / 2, cy / 2 + cyFace * 0.007, cxMouth, cyMouth, 1, (180 * 16) - 1);
			// painter.drawArc(1, 1, 100, 100, 1, (180 * 16) - 2);
		}
		else
		{
			painter.drawArc((cx - cxMouth) / 2, cy / 2 + cyFace * 0.007, cxMouth, cyMouth, (180 * 16) + 2, (180 * 16) - 2);
			// painter.drawArc(1, 1, 100, 100, (180 * 16) + 2, (180 * 16) - 2);
		};
#endif

	}


	// Now test the spectrum:
	for (int i88 = -1; i88 < 12; i88++) {
		double f;
		f = 0.1 * (double)i88;
		CRGBA mujKolor = K3_color(f, 1.0);
		int k3R, k3G, k3B, k3A;
		k3R = mujKolor.red(); k3G = mujKolor.green();
		k3B = mujKolor.blue(); k3A = mujKolor.alpha();
		QmujKolor.setAlpha((int)(mujKolor.alpha()));
		QmujKolor.setRed((int)(mujKolor.red()));
		QmujKolor.setGreen((int)(mujKolor.green()));
		QmujKolor.setBlue((int)(mujKolor.blue()));
		pen.setColor(QmujKolor);
		// painter.fillRect(10+i88*10, 2, 9, 30, QmujKolor);

	}

	// K3Face QQQ:
	QSize* K3Size = new(QSize);
	int K3cx, K3cy;
	*K3Size = this->size();
	K3cx = K3Size->width();
	K3cy = K3Size->height();

	// pen.setColor(QColor(255, 255, 0, 0));
	// painter.drawArc(QRect(2, 2, K3cx - 4, K3cy - 4), 16 * 30, 16 * 270);
	painter.end();


};

