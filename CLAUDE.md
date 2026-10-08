# pyTMDP - Claude 作業ルール

diphoton 質量スペクトルから top 質量を決める現象論コード（Python）と、そのテンプレートを作る Fortran（gg2aa, arXiv:1607.00990）。
物理の文脈は研究リポ `hiroshiyokoya/ai-coresearch`（スキル `hep-threshold-matching`）。**このリポは public。**

## 実装前に構想を確認すること

ユーザーの発言が提案・ディスカッションなのか、実装依頼なのかを判断する。曖昧なら実装せず方針を確認する。
実装・PR の作成は「やって」「実装して」「作って」など明示的に依頼されたときのみ。
`docs/REVIEW.md` の修正提案は、承認されたものだけを実装する。

## 正確性

- 数値・物理の主張には確度を付ける：**verified**（実行して確認）／**read**（コードを読んだ）／**inference**（推論）。
- 実行していないものを「動作確認済み」と言わない。テンプレートや断面積を変える変更は、変更前後の数値を並べて示す。
- 物理の仕様（カット・PDF・スケール・Green 関数）を変えるときは、コード上の修正と分けて議論する。

---

## イシュー管理ルール

### イシューを立てるときは必ずラベルを付ける

| ラベル | 用途 |
|--------|------|
| `bug` | 不具合修正 |
| `feature` | 新機能 |
| `enhancement` | 既存機能の改善・ドキュメント |
| `refactor` | リファクタリング |
| `infra` | Docker・ビルド環境 |
| `investigation` | 調査・レビュー |
| `maintenance` | 保守 |
| `tracking` | トラッキングイシュー |

マイルストーンも必ず付ける（`gh issue create --milestone`）。

### トラッキングイシューを常に最新に保つ

`tracking` ラベルのイシュー（現在: **#7**）は Claude が常に最新に保つ。マイルストーンごとにセクションを分け、
完了したイシューも履歴として残す。本文にアイコンやチェックボックスは付けない。

### 🚨 `gh issue create` とトラッキング更新は**セットの操作**

```bash
gh issue create -R hiroshiyokoya/pyTMDP --title "..." --label "..." --milestone "..." --body "..."
bash "$HOME/.claude/skills/git-dev-setup/scripts/update_tracking.sh" hiroshiyokoya/pyTMDP 7 "<見出し>" "- #<番号> タイトル"
```

---

## ブランチ・PR のルール

- **既定ブランチ = `develop`（開発本流）／ `master` = リリースブランチ**（2018 年版はここ）
- 通常開発: `issue` → `feature/<番号>-<簡潔な名前>` → `PR` → `develop`
- `develop` → `master` はユーザーが手動でマージ
- PR タイトルにイシュー番号（例: `feat: xxx (#12)`）、本文に `Closes #<番号>`
- **PR のマージはユーザーが行う。** Claude は `gh pr merge` を実行しない
- コミット前に `git -C D:/Physics/pyTMDP branch --show-current` でブランチ確認

---

## 環境

- すべて Docker イメージ `pytmdp:dev`（`docker/Dockerfile`）の中で動かす。リポを `/work` にマウントする。
  ```bash
  docker build -t pytmdp:dev docker
  MSYS_NO_PATHCONV=1 docker run --rm -v "D:/Physics/pyTMDP:/work" pytmdp:dev make -C fortran
  ```
  Git Bash から呼ぶときは `MSYS_NO_PATHCONV=1` を付ける（`/work` がパス変換されるのを防ぐ）。
- テスト：`python3 -m pytest -q tests`（コンテナ内）。
- 本番テンプレート 1 枚は 1 コアで約 80 分（低統計の実測からの見積もり）。長い生成ジョブは完了枚数を数えて見張る。
