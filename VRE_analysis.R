# =========================================================================
# A병원 VRE(Vancomycin-Resistant Enterococci) 감염 추세 분석
# 분석 모델: Poisson-Gamma Conjugate Model & Monte Carlo Simulation
# =========================================================================

# 1. 데이터 입력 (보정된 수치 사용)
# n: 총 재원일수 / 1000 (단위 노출량)
# y: 실제 발생 건수 (Rate * n)
n_2023 <- 261.2; y_2023 <- 5.2   # 2023년 데이터
n_2025 <- 334.0; y_2025 <- 110.2 # 2025년 데이터

# 2. 베이지안 사후분포 도출 (Conjugate Prior: Gamma Distribution)
# 무정보적 사전분포(Non-informative Prior) Gamma(0.001, 0.001) 가정
# Posterior ~ Gamma(alpha + y, beta + n)
alpha_prior <- 0.001
beta_prior <- 0.001

# 사후분포 파라미터 계산
post_shape_23 <- alpha_prior + y_2023
post_rate_23  <- beta_prior + n_2023

post_shape_25 <- alpha_prior + y_2025
post_rate_25  <- beta_prior + n_2025

# 3. 몬테카를로 시뮬레이션 (Monte Carlo Simulation)
# 사후분포에서 각각 100,000개의 샘플을 추출하여 비교
set.seed(42) # 재현성을 위한 시드 설정
samples_2023 <- rgamma(100000, shape = post_shape_23, rate = post_rate_23)
samples_2025 <- rgamma(100000, shape = post_shape_25, rate = post_rate_25)

# 2025년 발생률이 2023년보다 높을 확률 계산
prob_increase <- mean(samples_2025 > samples_2023)

# 4. 결과 시각화
x_range <- seq(0, 0.5, length.out = 1000)
dist_2023 <- dgamma(x_range, shape = post_shape_23, rate = post_rate_23)
dist_2025 <- dgamma(x_range, shape = post_shape_25, rate = post_rate_25)

plot(x_range, dist_2025, type = "l", col = "#E46432", lwd = 3,
     main = "VRE Infection Rate: Posterior Distribution Comparison",
     xlab = "Infection Rate (per 1,000 Patient-days)",
     ylab = "Density", xlim = c(0, 0.5))
lines(x_range, dist_2023, col = "#3B7EB9", lwd = 3)
abline(v = mean(samples_2023), col = "#3B7EB9", lty = 2)
abline(v = mean(samples_2025), col = "#E46432", lty = 2)
legend("topright", 
       legend = c(paste0("2023 (Mean: ", round(mean(samples_2023), 3), ")"),
                  paste0("2025 (Mean: ", round(mean(samples_2025), 3), ")")),
       col = c("#3B7EB9", "#E46432"), lwd = 3)

# 5. 수치 결과 출력
cat("--- 분석 결과 ---\n")
cat("2023년 95% 신용구간:", quantile(samples_2023, c(0.025, 0.975)), "\n")
cat("2025년 95% 신용구간:", quantile(samples_2025, c(0.025, 0.975)), "\n")
cat("2025년 발생률이 더 높을 사후 확률:", prob_increase * 100, "%\n")