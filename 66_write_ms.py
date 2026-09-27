import sys
content = open(r"D:\BCBM_Project\code\manuscript_content.txt", encoding="utf-8").read()
with open(r"D:\BCBM_Project\manuscript\manuscript_v1.md", "w", encoding="utf-8") as f:
    f.write(content)
print(f"Written {len(content)} chars")
