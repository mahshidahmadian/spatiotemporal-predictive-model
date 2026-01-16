---

## 📊 Data Schema & Relationships

### 1. Fish Detection Data (`fish_data.csv`)

| Column | Data Type | Description |
| :--- | :--- | :--- |
| `tagname` | `string` | Unique alphanumeric identifier for each fish. |
| `fish_id` | `integer` | Normalized numeric ID for computational efficiency. |
| `receiver_id` | `int / NA` | ID of the station that recorded the ping. `NA` indicates a missing timestep. |
| `longitude` | `num / NA` | Geographic longitude of the detection (Decimal Degrees). |
| `latitude` | `num / NA` | Geographic latitude of the detection (Decimal Degrees). |
| `date` | `datetime` | Timestamp of the detection event. |
| `miss_id` | `binary` | **Indicator:** `0` = Observed detection; `1` = Missing/Hidden state. |

### 2. Receiver Location Data (`receiver_data.csv`)

| Column | Data Type | Description |
| :--- | :--- | :--- |
| `receiver_id` | `integer` | Unique identifier matching the `fish_data` file. |
| `longitude` | `numeric` | Static longitude of the physical receiver station. |
| `latitude` | `numeric` | Static latitude of the physical receiver station. |

---
