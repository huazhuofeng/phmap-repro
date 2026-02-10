# Example

这个目录用于存放本项目的可复用示例。

## 文件说明

- `config.example.yaml`：示例配置文件（从主配置模板复制）

## 如何使用示例配置

1. 复制并修改输入路径：

```bash
cp example/config.example.yaml config/config.yaml
# 编辑 config/config.yaml 中 input/outdir 等参数
```

2. 运行 Snakemake：

```bash
snakemake --use-conda -j 4
```

## 建议

- 每次新任务建议复制一份配置为新名字（如 `config_run_lung.yaml`）并记录参数。
- `outdir` 建议按日期或批次命名，便于追踪不同运行。
