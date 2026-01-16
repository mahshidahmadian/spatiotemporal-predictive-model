---

## 📊 Data Format

### Fish Detection Data (`fish_data.csv`)

| Column | Type | Description |
|--------|------|-------------|
| `tagname` | string | Unique fish identifier |
| `fish_id` | integer | Numeric fish ID |
| `receiver_id` | integer/NA | Receiver where detected (NA if missing) |
| `longitude` | numeric/NA | Detection longitude (NA if missing) |
| `latitude` | numeric/NA | Detection latitude (NA if missing) |
| `date` | date | Detection date |
| `miss_id` | 0/1 | Binary: 0=observed, 1=missing |

### Receiver Location Data (`receiver_data.csv`)

| Column | Type | Description |
|--------|------|-------------|
| `receiver_id` | integer | Unique receiver identifier |
| `longitude` | numeric | Receiver longitude |
| `latitude` | numeric | Receiver latitude |

---
