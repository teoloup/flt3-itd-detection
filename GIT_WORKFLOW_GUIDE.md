# 🚀 Git Workflow Cheat Sheet for Beginners

## The Golden Rule: Always check status first!
```bash
git status
```

## 📋 Daily Workflow (4 Simple Steps)

### 1️⃣ CHANGE
- Edit your Python files in VS Code
- Save the files (Ctrl+S)

### 2️⃣ CHECK & ADD
```bash
# See what changed
git status

# Add specific files
git add filename.py

# OR add all changed files
git add .
```

### 3️⃣ COMMIT
```bash
git commit -m "your message here"
```

### 4️⃣ PUSH
```bash
git push
```

## 💡 Commit Message Examples

### Good commit messages:
```bash
git commit -m "fix: correct CIGAR parsing logic"
git commit -m "feat: add new validation method"
git commit -m "docs: update README installation steps"
git commit -m "test: add unit tests for consensus generation"
git commit -m "refactor: improve performance in ITD detection"
```

### Message prefixes:
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation
- `test:` - Adding tests
- `refactor:` - Code improvement
- `perf:` - Performance improvement

## 🔍 Useful Commands

### Check what's different
```bash
git diff filename.py  # See exact changes in a file
git log --oneline     # See recent commits
```

### Undo changes (BEFORE committing)
```bash
git checkout -- filename.py  # Undo changes to a file
git reset                     # Unstage all files
```

### Pull latest changes from GitHub
```bash
git pull  # Get updates from GitHub (do this before starting work)
```

## 🛡️ Safety Tips

1. **Always `git status` first** - know what you're doing
2. **Commit often** - small, frequent commits are better
3. **Pull before you push** - get latest changes first
4. **Write clear messages** - your future self will thank you
5. **Test before committing** - make sure your code works

## 🚨 Emergency Commands

### Oops, I committed something wrong!
```bash
git reset --soft HEAD~1  # Undo last commit, keep changes
```

### I want to see what changed
```bash
git show  # Show last commit details
```

### I want to go back to a previous version
```bash
git log --oneline        # Find the commit ID
git checkout COMMIT_ID   # Go to that version (temporary)
git checkout main        # Come back to latest
```

## 📱 Quick Reference Card

| What I want to do | Command |
|-------------------|---------|
| See what changed | `git status` |
| Add all changes | `git add .` |
| Add one file | `git add filename.py` |
| Commit changes | `git commit -m "message"` |
| Send to GitHub | `git push` |
| Get from GitHub | `git pull` |
| See commit history | `git log --oneline` |
| Undo last commit | `git reset --soft HEAD~1` |

## 🎯 Real Example Workflow

```bash
# 1. Start your day
git status
git pull

# 2. Make changes to your code
# ... edit files in VS Code ...

# 3. Check and commit
git status
git add .
git commit -m "fix: improve ITD validation accuracy"
git push

# 4. Done! ✅
```

## 🏠 Moving Your Project & Repository Connection

### How Git Knows Which Repository to Use
```bash
git remote -v  # Shows which GitHub repo you're connected to
```
Output: `origin https://github.com/teoloup/flt3-itd-detection`

### Moving Your Project Folder
✅ **SAFE**: Move entire project folder anywhere
- The hidden `.git` folder moves with it
- All Git history and connections preserved
- Continue using same commands

❌ **UNSAFE**: Copy only Python files without `.git` folder
- Loses all Git history and GitHub connection
- Would need to re-initialize Git

### If You Want to Change GitHub Repository
```bash
# See current connection
git remote -v

# Change to different repository
git remote set-url origin https://github.com/username/new-repo.git

# Verify the change
git remote -v
```

### Working from Multiple Computers
```bash
# On new computer, clone your repository
git clone https://github.com/teoloup/flt3-itd-detection.git
cd flt3-itd-detection

# Make changes normally
git add .
git commit -m "changes from laptop"
git push

# On original computer, get the changes
git pull
```

---
**Remember**: Practice makes perfect! Start with small changes and you'll get comfortable quickly. 🚀
