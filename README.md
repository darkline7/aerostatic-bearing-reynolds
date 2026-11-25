# Giải phương trình Reynolds 1D cho ổ trục khí tĩnh (Aerostatic Bearing)

Mô phỏng lại mô hình 1D của Miettinen et al. (2024), Precision Engineering 89 (2024) 1–10,
mã bài tập lớn môn Phương Pháp Số (ME2091), ĐH Bách Khoa TP.HCM.

## Mô hình

Phương trình Reynolds biến đổi trong tọa độ trụ (phương trình (6) của bài báo):

```
(1/r) d/dr( p r h^3 dp/dr ) = (6 k / h_p) (p^2 - p_s^2)
```

Điều kiện biên:

- Tại tâm r = 0: dp/dr = 0 (đối xứng). Về số, bắt đầu từ r_min = 1 µm để tránh điểm kỳ dị 1/r.
- Tại mép ngoài r = R = D/2 = 18.415 mm: p = p_a = 0.1 MPa.

Đưa về hệ phương trình vi phân bậc 1 với y = [p, dp/dr]:

```
dy0/dr = y1
dy1/dr = (6 k / (h^3 h_p y0)) (y0^2 - p_s^2) - (1/r + y1/y0) y1
```

Giải bằng `scipy.integrate.solve_bvp` (phương pháp collocation) trên 4 lưới
N = 10, 100, 1000, 10000 điểm, với 3 mức áp cấp p_s = 0.6, 0.4, 0.2 MPa.

## Tham số (bảng của đề bài)

| Tham số | Giá trị |
| --- | --- |
| Chiều dày màng khí h | 10 µm |
| Bán kính ổ R = D1/2 | 18.415 mm |
| Chiều dày vật liệu xốp h_p | 4.5 mm |
| Áp suất môi trường p_a | 0.1 MPa |
| Độ nhớt động lực μ | 1.78e-5 N·s/m² |
| Độ thấm k | 1.44e-15 m² (có thể chỉnh về 1.20e-15 bằng `--k`) |
| Áp suất cấp p_s | 0.6; 0.4; 0.2 MPa |

## Cài đặt và chạy

```bash
pip install -r requirements.txt
python3 solve_reynolds.py                 # mặc định k = 1.44e-15 m^2
python3 solve_reynolds.py --k 1.20e-15    # giá trị k đã hiệu chỉnh theo thực nghiệm
python3 solve_reynolds.py --paper-csv du_lieu_bai_bao.csv
```

File `du_lieu_bai_bao.csv` (tùy chọn) là dữ liệu số hóa từ đồ thị bài báo với 3 cột
`r_mm,p_MPa,ps_MPa` (ví dụ bằng WebPlotDigitizer), dùng để vẽ đường tam giác so sánh
trên hình phân bố áp suất.

## Kết quả đầu ra (thư mục `outputs/`)

- `fig1_pressure_distribution.png`: phân bố áp suất p(r) theo 3 mức p_s (N = 10000).
- `fig2_mesh_convergence.png`: so sánh 4 lưới N = 10/100/1000/10000 cho từng p_s.
- `fig3_error_vs_mesh.png`: sai số lớn nhất so với N = 10000 theo số điểm lưới.
- `pressure_Ps_*MPa_N10000.csv`: dữ liệu p(r) ở N = 10000 cho từng mức áp cấp.

Bảng in ra màn hình gồm áp suất tại tâm p(0) và khả năng tải
W = 2π ∫ (p − p_a) r dr (N).
