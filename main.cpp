#include <bits/stdc++.h>
using namespace std;

/*
 Needleman–Wunsch global alignment (simple linear gap penalty).
 - Inputs: two sequences (strings)
 - Parameters: match_score, mismatch_penalty (negative or 0), gap_penalty (negative)
 - Outputs: best alignment score and one optimal alignment (with '-' for gaps)
 Complexity: O(n * m) time and O(n * m) memory (can be optimized to linear memory
 for score only; but traceback requires the matrix).
*/

struct AlignmentResult {
    int score;
    string aligned_a;
    string aligned_b;
};

AlignmentResult needleman_wunsch(const string &a, const string &b,
                                 int match = 1, int mismatch = -1, int gap = -2,
                                 bool show_matrix = false) {
    int n = (int)a.size();
    int m = (int)b.size();
    // DP matrix: (n+1) x (m+1)
    vector<vector<int>> dp(n+1, vector<int>(m+1, INT_MIN));
    // direction matrix for traceback: 0 = diag, 1 = up (gap in b), 2 = left (gap in a)
    vector<vector<int>> dir(n+1, vector<int>(m+1, -1));

    // initialize
    dp[0][0] = 0;
    dir[0][0] = -1;
    for (int i = 1; i <= n; ++i) {
        dp[i][0] = dp[i-1][0] + gap;
        dir[i][0] = 1; // came from up (i-1,0)
    }
    for (int j = 1; j <= m; ++j) {
        dp[0][j] = dp[0][j-1] + gap;
        dir[0][j] = 2; // came from left (0,j-1)
    }

    // fill
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int score_diag = dp[i-1][j-1] + (a[i-1] == b[j-1] ? match : mismatch);
            int score_up   = dp[i-1][j] + gap;   // gap in b
            int score_left = dp[i][j-1] + gap;   // gap in a

            int best = score_diag;
            int best_dir = 0;
            if (score_up > best) {
                best = score_up;
                best_dir = 1;
            }
            if (score_left > best) {
                best = score_left;
                best_dir = 2;
            }
            dp[i][j] = best;
            dir[i][j] = best_dir;
        }
    }

    if (show_matrix) {
        cout << "DP matrix (scores):\n    ";
        for (int j = 0; j <= m; ++j) {
            if (j == 0) cout << "- ";
            else cout << b[j-1] << ' ';
        }
        cout << '\n';
        for (int i = 0; i <= n; ++i) {
            if (i == 0) cout << "- ";
            else cout << a[i-1] << ' ';
            for (int j = 0; j <= m; ++j) {
                cout << setw(3) << dp[i][j] << ' ';
            }
            cout << '\n';
        }
        cout << '\n';
    }

    // traceback from dp[n][m]
    int i = n, j = m;
    string aligned_a_rev, aligned_b_rev;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dir[i][j] == 0) {
            aligned_a_rev.push_back(a[i-1]);
            aligned_b_rev.push_back(b[j-1]);
            --i; --j;
        } else if (i > 0 && dir[i][j] == 1) {
            aligned_a_rev.push_back(a[i-1]);
            aligned_b_rev.push_back('-');
            --i;
        } else if (j > 0 && dir[i][j] == 2) {
            aligned_a_rev.push_back('-');
            aligned_b_rev.push_back(b[j-1]);
            --j;
        } else {
            // fallback (shouldn't normally occur): prefer diag then up then left
            if (i > 0 && j > 0) {
                aligned_a_rev.push_back(a[i-1]);
                aligned_b_rev.push_back(b[j-1]);
                --i; --j;
            } else if (i > 0) {
                aligned_a_rev.push_back(a[i-1]);
                aligned_b_rev.push_back('-');
                --i;
            } else {
                aligned_a_rev.push_back('-');
                aligned_b_rev.push_back(b[j-1]);
                --j;
            }
        }
    }
    reverse(aligned_a_rev.begin(), aligned_a_rev.end());
    reverse(aligned_b_rev.begin(), aligned_b_rev.end());

    AlignmentResult res;
    res.score = dp[n][m];
    res.aligned_a = move(aligned_a_rev);
    res.aligned_b = move(aligned_b_rev);
    return res;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cout << "Needleman–Wunsch Global Alignment (linear gap penalty)\n";
    cout << "Enter sequence A: ";
    string a, b;
    if (!getline(cin, a)) return 0;
    cout << "Enter sequence B: ";
    if (!getline(cin, b)) return 0;

    // optional: allow user to set scoring parameters
    int match = 1, mismatch = -1, gap = -2;
    cout << "Use default scoring? match=+1 mismatch=-1 gap=-2 (y/n): ";
    string yn;
    getline(cin, yn);
    if (!yn.empty() && (yn[0] == 'n' || yn[0] == 'N')) {
        cout << "Enter match score (int): ";
        string s;
        getline(cin, s);
        if (!s.empty()) match = stoi(s);
        cout << "Enter mismatch penalty (int, typically negative or 0): ";
        getline(cin, s);
        if (!s.empty()) mismatch = stoi(s);
        cout << "Enter gap penalty (int, typically negative): ";
        getline(cin, s);
        if (!s.empty()) gap = stoi(s);
    }

    cout << "Show DP matrix? (y/n): ";
    getline(cin, yn);
    bool show_matrix = (!yn.empty() && (yn[0] == 'y' || yn[0] == 'Y'));

    AlignmentResult result = needleman_wunsch(a, b, match, mismatch, gap, show_matrix);

    cout << "\nAlignment score: " << result.score << "\n\n";
    cout << "Aligned sequences:\n";
    cout << result.aligned_a << '\n';
    cout << result.aligned_b << '\n';

    // optionally show a matching line with '|' for matches and ' ' for mismatches/gaps
    string mid;
    for (size_t k = 0; k < result.aligned_a.size(); ++k) {
        if (result.aligned_a[k] == result.aligned_b[k]) mid.push_back('|');
        else mid.push_back(' ');
    }
    cout << mid << '\n';

    return 0;
}
