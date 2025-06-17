#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <limits>

using namespace std;

struct TableEntry {
    double x;
    double value;
};

vector<TableEntry> table4;
vector<TableEntry> table5;
map<string, double> table6;

bool loadTable4(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) return false;

    table4.clear();
    string line;
    getline(file, line);
    getline(file, line);

    TableEntry entry;
    while (file >> entry.x >> entry.value) {
        table4.push_back(entry);
    }

    return true;
}

bool loadTable5(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) return false;

    table5.clear();
    string line;
    getline(file, line);
    getline(file, line);

    TableEntry entry;
    while (file >> entry.x >> entry.value) {
        table5.push_back(entry);
    }

    return true;
}

bool loadTable6(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) return false;

    table6.clear();
    string line;
    getline(file, line);
    getline(file, line);

    string text;
    double x;
    while (file >> text >> x) {
        table6[text] = x;
    }

    return true;
}

double Min(double a, double b, double c, double d) {
    return min({ a, b, c, d });
}

double U1(double x) {
    if (fabs(x) < 1e-10) return 1.0; // Заміна 0 на 1 для уникнення ділення на 0
    return atan(asin(sin(3 * x)));
}

double T1(double x) {
    return atan(acos(sin(2 * x)));
}

double Qqn1(double x, double y, double z) {
    double ux = U1(x);
    return x / ux + y * T1(y) - U1(z) * T1(z);
}

double Qnk1(double x, double y) {
    return 1.15 * Qqn1(x, y, x + y) - 0.95 * Qqn1(y, x, x - y);
}

double Qnk_algorithm2(double x, double y) {
    return x * Qnk1(x, y) + y * Qnk1(y, x) - 0.05 * Qnk1(x, y) * Qnk1(y, x);
}

double GetFrom(const string& text) {
    if (table6.empty()) return 0.0;
    auto it = table6.find(text);
    return (it != table6.end()) ? it->second : 0.0;
}

double Stext(double x, const string& text) {
    if (text.empty()) {
        return GetFrom("tet") + GetFrom("set") - x;
    }
    else if (x <= 0) {
        return GetFrom("get") + GetFrom(text);
    }
    else {
        return GetFrom(text) + x;
    }
}

double Ktext(double x, double y, double z, const string& text) {
    double min_val = (z < 0) ? Min(x, y, x - z, y - z) : Min(x, y, z - x, z - y);
    double base = GetFrom(text.empty() ? "set" : text);

    // Уточнений діапазон [0.3, 0.7] з плавнішою характеристикою
    return 0.5 + 0.2 * tanh(0.5 * min_val + 0.3 * base - 0.7);
}

double interpolate(const vector<TableEntry>& table, double x) {
    if (table.empty()) return 0.0;
    if (x <= table.front().x) return table.front().value;
    if (x >= table.back().x) return table.back().value;

    for (size_t i = 0; i < table.size() - 1; ++i) {
        if (table[i].x <= x && x <= table[i + 1].x) {
            double x1 = table[i].x, x2 = table[i + 1].x;
            double y1 = table[i].value, y2 = table[i + 1].value;
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }
    }
    return 0.0;
}

double U(double x) {
    if (table4.empty() || fabs(x) <= 5.0) return U1(x);
    return interpolate(table4, x);
}

double T(double x) {
    if (table5.empty() || fabs(x) <= 10.0) return T1(x);
    return interpolate(table5, x);
}

double Qkn(double x, double y) {
    double ux = U(x);
    if (fabs(ux) < 1e-10) return 0.0;
    return x / ux + y * T(y);
}

double Qnk(double x, double y) {
    x = max(min(x, 5.0), -5.0);  // Жорсткіше обмеження
    y = max(min(y, 10.0), -10.0);

    if (fabs(x) < 0.001 && fabs(y) < 0.001) return 1.0;  // Захист від нулів

    return Qkn(x, y) + x * Qkn(y, x);
}

double Rsv(double x, double y, double z) {
    // Обмеження вхідних значень
    x = max(min(x, 5.0), -5.0);
    y = max(min(y, 5.0), -5.0);
    z = max(min(z, 5.0), -5.0);

    // Спрощені формули з обмеженням результатів
    if (z > x && z > y) return max(min(z * 0.5 * Qnk(x, y) - x * y, 50.0), -50.0);
    if (x > y && x > z) return max(min(x * 0.5 * Qnk(z, y) + y * z, 50.0), -50.0);
    if (y > x && y > z) return max(min(y * 0.5 * Qnk(x, z) + x * z, 50.0), -50.0);
    return 0.0;
}

double func(double x, double y, double z) {
    double rsv1 = Rsv(x, y, z);
    double rsv2 = Rsv(y, z, x);
    double rsv3 = Rsv(z, x, y);

    // Додаткове обмеження
    return max(min(rsv1 + 0.05 * rsv2 * rsv3, 100.0), -100.0);
}

double Tsm(double x, double y) {
    x = max(min(x, 6.0), -6.0);
    y = max(min(y, 6.0), -6.0);

    double term1 = pow(x, 4) - 3 * pow(x, 2) + 2;
    double term2 = 4 * pow(y, 4) - pow(x, 2);

    if (term2 <= 0 || term1 <= 0) return 0.02;

    // Суперстабільне обчислення
    return 0.08 * log(1 + 0.4 * fabs(term1 * term2));
}

double Mts(double x, double y) {
    double term = 4 * pow(y, 4) - pow(x, 2);
    if (term < 0) return 0.0;

    double tsm1 = Tsm(x, y);
    double tsm2 = Tsm(x, x);

    // Обмеження результатів Tsm
    tsm1 = max(min(tsm1, 1e4), -1e4);
    tsm2 = max(min(tsm2, 1e4), -1e4);

    return x * tsm1 - y * tsm2;
}

double Mis(double x, double y) {
    return Mts(x, y);
}

double Mlt(double x, double y, double z) {
    double term1 = x * max(Mis(x, y), 0.01);  // Гарантуємо мінімум 0.01
    double term2 = z * max(Mis(z, y), 0.01);
    return max(term1 + term2, 0.01);  // Мінімальне m = 0.01
}

double Y(double x) {
    if (100 - x * x < 0) return 0.0;
    double term = x * sqrt(100 - x * x);
    if (term <= 0) return 0.0;
    return (term < 1) ? 1.0 : log(term);
}

double func_regr(double r, double k, double m) {
    // Жорсткі обмеження
    r = max(min(r, 50.0), -50.0);
    k = max(min(k, 0.5), 0.1);
    m = max(min(m, 0.5), 0.01);

    // Остаточна формула
    return (2.0 * k * r) - (0.2 * m * r);
}

int main() {
    if (!loadTable4("dat1.dat"))
        cout << "Warning: Using Algorithm 2 for U(x)" << endl;
    if (!loadTable5("dat2.dat"))
        cout << "Warning: Using Algorithm 2 for T(x)" << endl;
    if (!loadTable6("dat3.dat"))
        cout << "Warning: GetFrom will return 0" << endl;

    double x, y, z;
    string text;

    cout << "Enter x y z values: ";
    cin >> x >> y >> z;
    cin.ignore();
    cout << "Enter text (empty if none): ";
    getline(cin, text);

    try {
        double r = func(x, y, z);
        double k = Ktext(x, y, z, text);
        double m = Mlt(x, y, z);

        // Додаткові обмеження перед фінальним обчисленням
        r = max(min(r, 1e4), -1e4);
        k = max(min(k, 1e4), -1e4);
        m = max(min(m, 1e4), -1e4);

        cout << "\nIntermediate values:" << endl;
        cout << "r = " << r << endl;
        cout << "k = " << k << endl;
        cout << "m = " << m << endl;

        double result = func_regr(r, k, m);
        cout << "\nFinal result: " << result << endl;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}