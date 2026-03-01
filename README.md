# Asynchronous Communication Distributed Multichannel Active Noise Control (AC-DMCANC)

The project aims to enable **communication-efficient**, **robust**, and **scalable** distributed ANC suitable for real-world deployment with heterogeneous network conditions.

---

## 📌 Overview

Distributed Multichannel Active Noise Control (DMCANC) achieves global noise reduction by distributing computation across multiple low-cost nodes. However, most existing methods assume **synchronous communication at every sampling instant**, which is impractical due to bandwidth and hardware constraints.

This repository implements an **asynchronous communication framework** where nodes:

- Operate independently using local adaptation
- Communicate only when performance degrades
- Exchange compact weight-difference information
- Maintain stability during non-communication phases


## ✨ Key Features

- 🔁 **Asynchronous communication** to reduce communication overhead  
- 🧱 **Weight-Constrained FxLMS (WCFxLMS)** for stability without communication  
- 🔄 **Mixed Weight Difference (MWD)** for efficient global information fusion  
- 📡 Event-triggered communication based on residual noise level  
- 🌐 Robust operation in heterogeneous and bandwidth-limited networks  

---


---

## ⚙️ How It Works

### Non-communication phase
- Each node runs **WCFxLMS**
- Filter updates are constrained to prevent divergence
- Local noise reduction is maintained

### Communication trigger
A node requests communication when:
- Its averaged residual noise level stops improving

### Communication phase
- Nodes exchange **weight differences**
- MWD combines global information
- Control filters and center points are updated

---

## 📊 Advantages

- ✔ Reduced communication frequency  
- ✔ Stable operation without continuous data exchange  
- ✔ Suitable for real-world hardware and wireless networks  

---

## ⚠️ Trade-offs

- Slightly slower convergence than fully synchronous methods  
- Slightly lower final noise reduction performance  
- Designed to balance performance and communication efficiency  

---

## Paper

https://arxiv.org/pdf/2601.15653
