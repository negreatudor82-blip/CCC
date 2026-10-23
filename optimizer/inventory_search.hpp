#ifndef MOLECULE_FACTORY_INVENTORY_SEARCH_HPP
#define MOLECULE_FACTORY_INVENTORY_SEARCH_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

// Exact uniform-cost search for small networks. Purchases are made only to
// cover the gross inputs of a reaction. Inventories are never clipped: hitting
// either resource limit returns the incumbent without claiming optimality.
namespace inventory_search {

using ll = long long;

struct Reaction {
    std::vector<std::pair<int, ll>> inputs, outputs;
};

struct Result {
    ll cost;
    std::vector<int> sequence;
    bool exhausted;
};

namespace detail {
struct InventoryHash {
    std::size_t operator()(const std::vector<ll>& inventory) const noexcept {
        std::size_t h = inventory.size();
        for (ll amount : inventory) {
            const std::size_t x = std::hash<ll>{}(amount);
            h ^= x + std::size_t(0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
        }
        return h;
    }
};

struct Node {
    std::vector<ll> inventory;
    ll cost;
    std::size_t parent;
    int reaction;
};

struct CompactReaction {
    int original;
    std::vector<std::pair<int, ll>> inputs, outputs;
};
}  // namespace detail

inline Result search(
    const std::vector<ll>& prices,
    const std::vector<Reaction>& reactions,
    int target,
    ll needed,
    ll incumbentCost,
    const std::vector<int>& incumbentSequence,
    std::chrono::steady_clock::time_point deadline,
    std::size_t maxStates = 100000) {

    Result answer{incumbentCost, incumbentSequence, false};
    const std::size_t noParent = std::numeric_limits<std::size_t>::max();
    const ll limit = std::numeric_limits<ll>::max();
    const int moleculeCount = static_cast<int>(prices.size());
    if (target < 0 || target >= moleculeCount) return answer;
    if (needed <= 0) return Result{0, {}, true};
    if (incumbentCost == 0 && !incumbentSequence.empty()) {
        answer.exhausted = true;
        return answer;
    }
    if (maxStates == 0 || std::chrono::steady_clock::now() >= deadline)
        return answer;

    // Backward closure is safe even with catalysts and cycles: whenever a
    // reaction may produce a required molecule, retain it and all its gross
    // inputs. Outputs outside this closure cannot enable any retained reaction.
    std::vector<std::vector<int>> producers(moleculeCount);
    for (int r = 0; r < static_cast<int>(reactions.size()); ++r)
        for (const auto& output : reactions[r].outputs)
            if (output.first >= 0 && output.first < moleculeCount && output.second > 0)
                producers[output.first].push_back(r);

    std::vector<unsigned char> relevantMolecule(moleculeCount, 0);
    std::vector<unsigned char> relevantReaction(reactions.size(), 0);
    std::vector<int> pending{target};
    relevantMolecule[target] = 1;
    for (std::size_t p = 0; p < pending.size(); ++p) {
        for (int r : producers[pending[p]]) {
            if (relevantReaction[r]) continue;
            relevantReaction[r] = 1;
            for (const auto& input : reactions[r].inputs) {
                const int m = input.first;
                if (m >= 0 && m < moleculeCount && input.second > 0 &&
                    !relevantMolecule[m]) {
                    relevantMolecule[m] = 1;
                    pending.push_back(m);
                }
            }
        }
    }
    std::vector<int> compactId(moleculeCount, -1);
    std::vector<ll> compactPrices;
    for (int m = 0; m < moleculeCount; ++m)
        if (relevantMolecule[m]) {
            compactId[m] = static_cast<int>(compactPrices.size());
            compactPrices.push_back(prices[m]);
        }
    const int compactTarget = compactId[target];
    const int activeCount = static_cast<int>(compactPrices.size());
    bool arithmeticLimited = false;
    std::vector<detail::CompactReaction> activeReactions;
    for (int r = 0; r < static_cast<int>(reactions.size()); ++r) {
        if (!relevantReaction[r]) continue;
        detail::CompactReaction reaction;
        reaction.original = r;
        bool valid = true;
        auto normalize = [&](const std::vector<std::pair<int, ll>>& source,
                             std::vector<std::pair<int, ll>>& destination,
                             bool isInput) {
            std::vector<__int128> totals(activeCount, 0);
            for (const auto& item : source) {
                if (item.first < 0 || item.first >= moleculeCount || item.second < 0) {
                    valid = false;
                    return;
                }
                int m = compactId[item.first];
                if (m < 0) {
                    if (isInput && item.second > 0) valid = false;
                    continue;
                }
                totals[m] += item.second;
                if (totals[m] > limit) {
                    valid = false;
                    arithmeticLimited = true;
                    return;
                }
            }
            for (int m = 0; m < activeCount; ++m)
                if (totals[m] > 0)
                    destination.emplace_back(m, static_cast<ll>(totals[m]));
        };
        normalize(reactions[r].inputs, reaction.inputs, true);
        normalize(reactions[r].outputs, reaction.outputs, false);
        if (valid) activeReactions.push_back(std::move(reaction));
    }

    std::vector<detail::Node> nodes;
    nodes.reserve(std::min<std::size_t>(maxStates, 100000));
    std::unordered_map<std::vector<ll>, std::size_t, detail::InventoryHash> seen;
    seen.reserve(std::min<std::size_t>(maxStates, 100000));
    std::vector<ll> empty(activeCount, 0);
    nodes.push_back(detail::Node{empty, 0, noParent, -1});
    seen.emplace(empty, 0);
    using QueueEntry = std::pair<ll, std::size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    queue.emplace(0, 0);
    auto reconstruct = [&](std::size_t id, int lastReaction) {
        std::vector<int> result;
        if (lastReaction >= 0) result.push_back(lastReaction);
        while (id != noParent && nodes[id].reaction >= 0) {
            result.push_back(nodes[id].reaction);
            id = nodes[id].parent;
        }
        std::reverse(result.begin(), result.end());
        return result;
    };

    while (!queue.empty()) {
        if (std::chrono::steady_clock::now() >= deadline) return answer;
        const QueueEntry entry = queue.top();
        queue.pop();
        const ll currentCost = entry.first;
        const std::size_t currentId = entry.second;
        if (currentCost != nodes[currentId].cost) continue;
        if (currentCost >= answer.cost) continue;

        // Copy before insertion: nodes may reallocate while children are added.
        const std::vector<ll> inventory = nodes[currentId].inventory;
        for (const auto& reaction : activeReactions) {
            if (std::chrono::steady_clock::now() >= deadline) return answer;
            __int128 nextCost = currentCost;
            bool legal = true;
            for (const auto& input : reaction.inputs) {
                const int m = input.first;
                const ll deficit = input.second > inventory[m] ? input.second - inventory[m] : 0;
                if (deficit == 0) continue;
                if (compactPrices[m] < 0) { legal = false; break; }
                nextCost += static_cast<__int128>(deficit) * compactPrices[m];
                if (nextCost >= answer.cost) { legal = false; break; }
            }
            if (!legal) continue;
            if (nextCost > limit) { arithmeticLimited = true; continue; }
            std::vector<ll> nextInventory = inventory;
            for (const auto& input : reaction.inputs)
                nextInventory[input.first] = std::max<ll>(0, nextInventory[input.first] - input.second);
            for (const auto& output : reaction.outputs) {
                const __int128 amount = static_cast<__int128>(nextInventory[output.first]) + output.second;
                if (amount > limit) { legal = false; arithmeticLimited = true; break; }
                nextInventory[output.first] = static_cast<ll>(amount);
            }
            if (!legal) continue;
            const ll candidateCost = static_cast<ll>(nextCost);
            if (nextInventory[compactTarget] >= needed) {
                answer.cost = candidateCost;
                answer.sequence = reconstruct(currentId, reaction.original);
                if (candidateCost == 0) { answer.exhausted = true; return answer; }
                continue;
            }
            auto old = seen.find(nextInventory);
            if (old != seen.end()) {
                const std::size_t id = old->second;
                if (candidateCost >= nodes[id].cost) continue;
                nodes[id].cost = candidateCost;
                nodes[id].parent = currentId;
                nodes[id].reaction = reaction.original;
                queue.emplace(candidateCost, id);
                continue;
            }
            if (nodes.size() >= maxStates) return answer;
            const std::size_t id = nodes.size();
            seen.emplace(nextInventory, id);
            nodes.push_back(detail::Node{std::move(nextInventory), candidateCost,
                                         currentId, reaction.original});
            queue.emplace(candidateCost, id);
        }
    }
    answer.exhausted = !arithmeticLimited;
    return answer;
}

}  // namespace inventory_search

#endif  // MOLECULE_FACTORY_INVENTORY_SEARCH_HPP
