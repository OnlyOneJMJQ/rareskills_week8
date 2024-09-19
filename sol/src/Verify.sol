// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.13;

contract Verify {
    uint256 constant PRIME_Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;

    struct G1Point {
        uint256 x;
        uint256 y;
    }

    struct G2Point {
        uint256[2] x;
        uint256[2] y;
    }
  
    function negate(G1Point memory p) internal pure returns (G1Point memory) {
        // The prime q in the base field F_q for G1
        if (p.x == 0 && p.y == 0) {
            return G1Point(0, 0);
        } else {
            return G1Point(p.x, PRIME_Q - (p.y % PRIME_Q));
        }
    }

    /// @dev Adds two EC points together and returns the resulting point.
    /// @param x1 The x coordinate of the first point
    /// @param y1 The y coordinate of the first point
    /// @param x2 The x coordinate of the second point
    /// @param y2 The y coordinate of the second point
    /// @return x The x coordinate of the resulting point
    /// @return y The y coordinate of the resulting point
    function ec_add(uint256 x1, uint256 y1, uint256 x2, uint256 y2) public view returns (uint256 x, uint256 y) {
        (bool ok, bytes memory result) = address(6).staticcall(abi.encode(x1, y1, x2, y2));
        require(ok, "add failed");
        (x, y) = abi.decode(result, (uint256, uint256));
    }

    /// @dev Multiplies an EC point by a scalar and returns the resulting point.
    /// @param scalar The scalar to multiply by
    /// @param x1 The x coordinate of the point
    /// @param y1 The y coordinate of the point
    /// @return x The x coordinate of the resulting point
    /// @return y The y coordinate of the resulting point
    function scalar_mul(uint256 scalar, uint256 x1, uint256 y1) internal view returns (uint256 x, uint256 y) {
        (bool ok, bytes memory result) = address(7).staticcall(abi.encode(x1, y1, scalar));
        require(ok, "mul failed");
        (x, y) = abi.decode(result, (uint256, uint256));
    }

    function pairing(
        G1Point memory a1,
        G2Point memory a2,
        G1Point memory b1,
        G2Point memory b2,
        G1Point memory c1,
        G2Point memory c2,
        G1Point memory d1,
        G2Point memory d2
    ) internal view returns (bool) {
        G1Point[4] memory p1 = [a1, b1, c1, d1];
        G2Point[4] memory p2 = [a2, b2, c2, d2];

        uint256 inputSize = 24;
        uint256[] memory input = new uint256[](inputSize);

        for (uint256 i = 0; i < 4; i++) {
            uint256 j = i * 6;
            input[j + 0] = p1[i].x;
            input[j + 1] = p1[i].y;
            input[j + 2] = p2[i].x[0];
            input[j + 3] = p2[i].x[1];
            input[j + 4] = p2[i].y[0];
            input[j + 5] = p2[i].y[1];
        }

        uint256[1] memory out;
        bool success;

        // solium-disable-next-line security/no-inline-assembly
        assembly {
            success := staticcall(sub(gas(), 2000), 8, add(input, 0x20), mload(inputSize), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }

        require(success, "pairing-opcode-failed");

        return out[0] != 0;
    }

    function verify(G1Point calldata A, G2Point calldata B, G1Point calldata Alpha, G2Point calldata Beta, G1Point[] calldata pTauPub, uint256[2] calldata pubInputs, G2Point calldata Gamma, G1Point calldata C, G2Point calldata Delta) public view returns (bool verified) {
        G1Point memory X;
        
        // Calculate x
        for (uint256 i = 0; i < pubInputs.length; i++) {
            if (pTauPub[i].x == 0 && pTauPub[i].y == 0) {
                continue;
            }

            if (i == 0) {
                (X.x, X.y) = scalar_mul(pubInputs[i], pTauPub[i].x, pTauPub[i].y);
            }
            else {
                G1Point memory temp;
                (temp.x, temp.y) = scalar_mul(pubInputs[i], pTauPub[i].x, pTauPub[i].y);
                (X.x, X.y) = ec_add(X.x, X.y, temp.x, temp.y);
            }
        }
        
        verified = pairing(negate(A), B, Alpha, Beta, X, Gamma, C, Delta);
    }
}
